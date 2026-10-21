import subprocess
import argparse
import sys
import re
import shutil
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent

STEPS = ['mafft', 'clean', 'qc_basic', 'pre_filter', 'iqtree1', 'raccoon', 'adar_mask', 'qc_phylo', 'iqtree2', 'temporal_root', 'filter_outliers', 'iqtree3', 'final_root']

REVIEW_POINT = 'temporal_root'


def run(cmd, description=None):
    if description:
        print(f"\n{'='*60}")
        print(f"  {description}")
        print(f"{'='*60}")
    print(f"$ {cmd}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"\nERROR: command failed with return code {result.returncode}")
        sys.exit(1)


def parse_best_model(iqtree_file):
    with open(iqtree_file) as f:
        for line in f:
            m = re.match(r'Best-fit model according to BIC:\s+(.+)', line)
            if m:
                return m.group(1).strip()
    print(f"ERROR: could not find best-fit model in {iqtree_file}")
    sys.exit(1)


def count_sequences(fasta_path):
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def filter_sequences(input_aln, output_aln, exclude_files):
    """Remove sequences whose ID or accession appears in any exclude file."""
    to_remove = set()
    for ef in exclude_files:
        if not Path(ef).exists():
            continue
        with open(ef) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                token = line.split()[0]
                if '|' in token:
                    to_remove.add(token.split('|')[0])
                    to_remove.add(token.split('|')[1].split('.')[0])
                else:
                    to_remove.add(token.split('.')[0])

    kept = 0
    removed = 0
    with open(output_aln, 'w') as out:
        with open(input_aln) as inp:
            writing = True
            for line in inp:
                if line.startswith('>'):
                    header = line[1:].strip()
                    fields = header.split('|')
                    seq_id = fields[0]
                    acc = fields[1].split('.')[0] if len(fields) > 1 else header.split('.')[0]
                    if seq_id in to_remove or acc in to_remove:
                        writing = False
                        removed += 1
                    else:
                        writing = True
                        kept += 1
                        out.write(line)
                elif writing:
                    out.write(line)
    print(f"Filtered: {kept} kept, {removed} removed")
    return kept


def find_phylo_flags(outdir):
    for f in sorted(Path(outdir).glob('*.phylo_flags.csv')):
        return str(f)
    print(f"ERROR: no phylo_flags.csv found in {outdir}")
    sys.exit(1)


def should_run(step, skip_to, stop_after):
    idx = STEPS.index(step)
    if skip_to is not None and idx < STEPS.index(skip_to):
        return False
    if stop_after is not None and idx > STEPS.index(stop_after):
        return False
    return True


parser = argparse.ArgumentParser(
    description="Run the Pathoplexus → ML tree pipeline",
    epilog="""Steps:
  [mafft] → clean → qc_basic → pre_filter → iqtree1 → raccoon → adar_mask → qc_phylo → iqtree2 → temporal_root
  [review outliers]
  filter_outliers → iqtree3 → final_root

Exclusion files:
  --exclude-file    Identity (non-consortium, same-patient duplicates): removed from everything
  --exclude-quality Quality (low-quality sequences): removed from phylogenetics only

By default, runs clean through temporal_root and stops for outlier review.
Use --skip-to to resume from any step. Use --stop-after to stop early.
Use --outlier-removal-threshold 3 to auto-remove outliers and run straight through.""",
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-f", "--fasta", required=True, help="FASTA file (pre-aligned from Pathoplexus, or unaligned if using --mafft)")
parser.add_argument("-m", "--metadata", required=True, help="Pathoplexus metadata TSV")
parser.add_argument("--mafft", action="store_true",
                    help="Run MAFFT alignment before cleaning (skip if sequences are already aligned from Pathoplexus)")
parser.add_argument("-p", "--prefix", required=True, help="Output prefix for cleaned files")
parser.add_argument("-a", "--aliases-file", default="data/aliases.csv",
                    help="Aliases CSV file")
parser.add_argument("--trim-start", type=int, default=None, help="First position to keep (1-based, inclusive)")
parser.add_argument("--trim-end", type=int, default=None, help="Last position to keep (1-based, inclusive)")
parser.add_argument("--skip-to", choices=STEPS, default=None,
                    help="Resume from this step (skips earlier steps, assumes their outputs exist)")
parser.add_argument("--stop-after", choices=STEPS, default=None,
                    help="Stop after this step (default: temporal_root unless --skip-to is past the review point)")
parser.add_argument("--adar-mode", choices=['mask', 'exclude_tc', 'exclude_all'], default='exclude_tc',
                    help="ADAR handling: mask sites, exclude T>C genomes, or exclude T>C/A>G/G>A genomes (default: exclude_tc)")
parser.add_argument("--seed", type=int, default=42,
                    help="Random seed for IQ-TREE reproducibility (default: 42)")
parser.add_argument("-d", "--outdir", default=None, help="Working directory for outputs (default: prefix parent directory)")
parser.add_argument("--exclude-file", type=str, default=None,
                    help="Identity exclusions (non-consortium, same-patient duplicates): excluded from everything")
parser.add_argument("--exclude-quality", type=str, default=None,
                    help="Quality exclusions (low-quality sequences): excluded from phylogenetics only")
parser.add_argument("--outlier-file", type=str, default=None,
                    help="Path to outlier exclusions file (default: {outdir}/outliers.txt)")
parser.add_argument("--outlier-removal-threshold", type=float, default=None,
                    help="Auto-remove outliers with |z| >= this threshold (e.g. 3.0) and skip manual review")
parser.add_argument("--model", type=str, default=None,
                    help="Use this substitution model (e.g. GTR+F+R3) instead of running ModelFinder")
parser.add_argument("--min-date", type=str, default=None,
                    help="Remove sequences before this date (YYYY-MM-DD), passed to clean_up_pathoplexus.py")
args = parser.parse_args()

outdir = Path(args.outdir) if args.outdir else Path(args.prefix).parent
outdir.mkdir(parents=True, exist_ok=True)

trim_args = ""
if args.trim_start is not None:
    trim_args += f" --trim_start {args.trim_start}"
if args.trim_end is not None:
    trim_args += f" --trim_end {args.trim_end}"

PYTHON = sys.executable
clean_script = SCRIPTS_DIR / "clean_up_pathoplexus.py"
mask_adar_script = SCRIPTS_DIR / "mask_adar_sites.py"
genome_qc_script = SCRIPTS_DIR / "genome_qc_report.py"
phylo_qc_script = SCRIPTS_DIR / "phylo_qc_report.py"
temporal_root_script = SCRIPTS_DIR / "temporal_root.py"
adar_file = outdir / "adar_exclusions.txt"
exclude_file = Path(args.exclude_file) if args.exclude_file else None
exclude_quality = Path(args.exclude_quality) if args.exclude_quality else None
outlier_file = Path(args.outlier_file) if args.outlier_file else outdir / "outliers.txt"
model_file = outdir / "best_model.txt"

# Determine stop point
if args.stop_after:
    stop_after = args.stop_after
elif args.outlier_removal_threshold is not None:
    stop_after = None  # auto-remove mode, run to end
elif args.skip_to and STEPS.index(args.skip_to) > STEPS.index(REVIEW_POINT):
    stop_after = None  # past review point, run to end
else:
    stop_after = REVIEW_POINT


def get_best_model():
    if model_file.exists():
        return open(model_file).read().strip()
    print(f"ERROR: {model_file} not found. Run iqtree1 step first.")
    sys.exit(1)


# --- Paths ---
adar_tag = 'ADAR_masked' if args.adar_mode == 'mask' else 'ADAR_excluded'

prefix_stem = Path(args.prefix).name
cleaned_prefix = outdir / f"{prefix_stem}.cleaned"
cleaned_aln = f"{cleaned_prefix}.aln"
base_prefix = outdir / f"{prefix_stem}"
base_aln = f"{base_prefix}.aln"
adar_prefix = outdir / f"{prefix_stem}.{adar_tag}"
adar_aln = f"{adar_prefix}.aln"
final_prefix = outdir / f"{prefix_stem}.{adar_tag}.final"
final_aln = f"{final_prefix}.aln"

# --- Steps ---

input_fasta = args.fasta
if should_run('mafft', args.skip_to, stop_after) and args.mafft:
    mafft_out = outdir / f"{prefix_stem}.mafft.fasta"
    run(f"mafft --auto --thread -1 {args.fasta} > {mafft_out}",
        "MAFFT alignment")
    input_fasta = str(mafft_out)
    print(f"MAFFT output: {count_sequences(input_fasta)} sequences")

if should_run('clean', args.skip_to, stop_after):
    min_date_arg = f" --min-date {args.min_date}" if args.min_date else ""
    run(f"{PYTHON} {clean_script} -f {input_fasta} -m {args.metadata} -a {args.aliases_file} -p {cleaned_prefix}{trim_args}{min_date_arg}",
        "Cleaning alignment (renaming, trimming, default exclusions)")
    print(f"Cleaned: {count_sequences(cleaned_aln)} sequences")

if should_run('qc_basic', args.skip_to, stop_after):
    qc_basic_report = outdir / "qc_report_basic.txt"
    run(f"{PYTHON} {genome_qc_script} -a {cleaned_aln} -m {args.metadata} -o {qc_basic_report}",
        "Basic genome QC report")

if should_run('pre_filter', args.skip_to, stop_after):
    exclude_lists = []
    if exclude_file and exclude_file.exists():
        exclude_lists.append(str(exclude_file))
    elif exclude_file:
        print(f"Warning: --exclude-file {exclude_file} not found")
    if exclude_quality and exclude_quality.exists():
        exclude_lists.append(str(exclude_quality))
    elif exclude_quality:
        print(f"Warning: --exclude-quality {exclude_quality} not found")

    if exclude_lists:
        n_before = count_sequences(cleaned_aln)
        n_after = filter_sequences(cleaned_aln, base_aln, exclude_lists)
        print(f"\n{'='*60}")
        print(f"  Pre-filter: {n_before} → {n_after} sequences ({n_before - n_after} excluded)")
        print(f"{'='*60}")
    else:
        shutil.copy2(cleaned_aln, base_aln)
        print(f"Base alignment: {count_sequences(base_aln)} sequences (no exclusions applied)")

if should_run('iqtree1', args.skip_to, stop_after):
    if args.model:
        run(f"iqtree3 -nt auto -s {base_aln} -m {args.model} -czb --ancestral -redo -seed {args.seed}",
            f"IQ-TREE with {args.model} + ancestral state reconstruction")
        with open(model_file, 'w') as f:
            f.write(args.model)
        print(f"\nModel: {args.model} (saved to {model_file})")
    else:
        run(f"iqtree3 -nt auto -s {base_aln} -m MFP -czb --ancestral -redo -seed {args.seed}",
            "IQ-TREE with ModelFinder + ancestral state reconstruction")
        best_model = parse_best_model(f"{base_aln}.iqtree")
        with open(model_file, 'w') as f:
            f.write(best_model)
        print(f"\nBest-fit model: {best_model} (saved to {model_file})")

treefile = f"{base_aln}.treefile"
statefile = f"{base_aln}.state"

if should_run('raccoon', args.skip_to, stop_after):
    run(f"raccoon tree-qc -t {treefile} --alignment {base_aln} --asr-state {statefile} --run-adar -d {outdir}",
        "Running raccoon tree-qc for ADAR identification")

if should_run('adar_mask', args.skip_to, stop_after):
    phylo_flags = find_phylo_flags(outdir)
    run(f"{PYTHON} {mask_adar_script} -a {base_aln} -r {phylo_flags} -t {treefile} -s {statefile} --mode {args.adar_mode} -p {adar_prefix}",
        f"ADAR handling (mode: {args.adar_mode})")
    print(f"After ADAR: {count_sequences(adar_aln)} sequences")

if should_run('qc_phylo', args.skip_to, stop_after):
    phylo_flags = find_phylo_flags(outdir)
    branch_snps = str(phylo_flags).replace('.phylo_flags.csv', '.branch_snps.reconstruction.csv')
    qc_basic_report = outdir / "qc_report_basic.txt"
    qc_phylo_report = outdir / "qc_report_phylo.txt"
    basic_qc_arg = f" --basic-qc {qc_basic_report}" if qc_basic_report.exists() else ""
    run(f"{PYTHON} {phylo_qc_script} -f {phylo_flags} -b {branch_snps} -m {args.metadata}{basic_qc_arg} -o {qc_phylo_report}",
        "Phylogenetic QC report")

if should_run('iqtree2', args.skip_to, stop_after):
    best_model = get_best_model()
    run(f"iqtree3 -nt auto -s {adar_aln} -m {best_model} -czb -redo -seed {args.seed}",
        f"IQ-TREE with {best_model}")

adar_tree = f"{adar_aln}.treefile"

if should_run('temporal_root', args.skip_to, stop_after):
    run(f"{PYTHON} {temporal_root_script} -t {adar_tree} --iterate --plot -o {outdir} --prefix {prefix_stem}.{adar_tag}",
        "Temporal rooting + outlier detection")
    detected_outlier_file = outdir / f"{prefix_stem}.{adar_tag}.outliers.txt"
    if detected_outlier_file.exists():
        if args.outlier_removal_threshold is not None:
            threshold = args.outlier_removal_threshold
            kept = []
            with open(detected_outlier_file) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    m = re.search(r'z=[+-]?([\d.]+)', line)
                    if m and float(m.group(1)) >= threshold:
                        kept.append(line)
            with open(outlier_file, 'w') as f:
                f.write(f"# Outliers auto-removed at |z| >= {threshold}\n")
                f.writelines(kept)
            print(f"Auto-removal: {len(kept)} outliers with |z| >= {threshold} (from {sum(1 for l in open(detected_outlier_file) if l.strip() and not l.startswith('#'))} detected at z>=2)")
        else:
            shutil.copy2(detected_outlier_file, outlier_file)
            n_outliers = sum(1 for line in open(outlier_file) if line.strip() and not line.startswith('#'))
            print(f"Outliers detected: {n_outliers} (copied to {outlier_file})")

if should_run('filter_outliers', args.skip_to, stop_after):
    if not Path(adar_aln).exists():
        print(f"ERROR: {adar_aln} not found. Run earlier steps first.")
        sys.exit(1)
    if not outlier_file.exists():
        print(f"ERROR: {outlier_file} not found. Run temporal_root step first, or create manually.")
        sys.exit(1)
    print(f"\n{'='*60}")
    print(f"  Filtering outliers")
    print(f"{'='*60}")
    n_final = filter_sequences(adar_aln, final_aln, [str(outlier_file)])
    print(f"Final: {n_final} sequences")

if should_run('iqtree3', args.skip_to, stop_after):
    best_model = get_best_model()
    run(f"iqtree3 -nt auto -s {final_aln} -m {best_model} -czb -abayes -redo -seed {args.seed}",
        f"Final IQ-TREE with {best_model} + aBayes support")

final_tree = f"{final_aln}.treefile"

if should_run('final_root', args.skip_to, stop_after):
    run(f"{PYTHON} {temporal_root_script} -t {final_tree} --plot -o {outdir} --prefix {prefix_stem}.{adar_tag}.final",
        "Final temporal rooting")

# --- Summary ---
last_step = stop_after or STEPS[-1]
if last_step == REVIEW_POINT:
    print(f"\n{'='*60}")
    print(f"  Stopped after temporal rooting. (ADAR mode: {args.adar_mode})")
    print(f"{'='*60}")
    print(f"\nRooted tree:       {outdir / f'{prefix_stem}.{adar_tag}.bestFitRoot.tree'}")
    print(f"RTT plot:          {outdir / f'{prefix_stem}.{adar_tag}.rtt_plot.pdf'}")
    print(f"Detected outliers: {outdir / f'{prefix_stem}.{adar_tag}.outliers.txt'} (original, do not edit)")
    print(f"Working copy:      {outlier_file} (edit this one)")
    print(f"\nReview outliers, edit {outlier_file} if needed, then continue:")
    print(f"  python {sys.argv[0]} ... --skip-to filter_outliers --adar-mode {args.adar_mode}")
elif STEPS.index(last_step) >= STEPS.index('final_root'):
    n_final = count_sequences(final_aln)
    n_outliers = sum(1 for line in open(outlier_file) if line.strip() and not line.startswith('#'))
    adar_info = ""
    if adar_file.exists():
        n_adar = sum(1 for line in open(adar_file) if line.strip() and not line.startswith('#'))
        adar_info = f", {n_adar} ADAR-excluded" if args.adar_mode.startswith('exclude') else ", ADAR sites masked"
    print(f"\n{'='*60}")
    print(f"  Pipeline complete.")
    print(f"{'='*60}")
    print(f"\nFinal alignment:   {final_aln}")
    print(f"Final rooted tree: {outdir / f'{prefix_stem}.{adar_tag}.final.bestFitRoot.tree'}")
    print(f"Final RTT plot:    {outdir / f'{prefix_stem}.{adar_tag}.final.rtt_plot.pdf'}")
    print(f"ADAR mode: {args.adar_mode}{adar_info}")
    print(f"Outlier exclusions: {n_outliers}")
    print(f"Final sequences: {n_final}")
