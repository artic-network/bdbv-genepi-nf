#!/usr/bin/env python3
"""End-to-end update: download → pipeline → health checks → BEAST XMLs.

Chains download_pathoplexus.py, run_pipeline.py, health_checks.py, and
generate_beast_xml.py into a single workflow.

Exclude lists (canonical, in data/):
  exclude_identity.txt  — confirmed non-consortium duplicates
  exclude_quality.txt   — garbage sequences
  exclude_pending.txt   — under review, temporarily excluded

Duplicate_sample flags from health checks are auto-appended to exclude_pending.txt.

Usage:
    python scripts/run_update.py \
        --previous-data data/20260915 \
        --template-sg templates/GTR_SG.template.xml \
        --template-egc templates/GTR_EGC.template.xml \
        --anchor-date 2025-11-08
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent
PYTHON = sys.executable


def read_exclude_entries(path):
    entries = set()
    if not Path(path).exists():
        return entries
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                entries.add(line.split()[0])
    return entries


def merge_exclude_files(output_path, *input_paths):
    """Merge multiple exclude files into one, deduplicating."""
    seen = set()
    lines = []
    for p in input_paths:
        if not Path(p).exists():
            continue
        with open(p) as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith('#'):
                    continue
                entry = stripped.split()[0]
                if entry not in seen:
                    seen.add(entry)
                    lines.append(entry)
    with open(output_path, 'w') as f:
        f.write('# Merged exclude list (identity + quality + pending)\n')
        for entry in lines:
            f.write(entry + '\n')
    return len(lines)


def run(cmd, description=None, capture=False):
    if description:
        print(f"\n{'='*60}")
        print(f"  {description}")
        print(f"{'='*60}")
    print(f"$ {cmd}\n")
    if capture:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(result.stdout, end='')
            print(result.stderr, file=sys.stderr, end='')
            print(f"\nERROR: command failed with return code {result.returncode}")
            sys.exit(1)
        print(result.stdout, end='')
        return result.stdout
    else:
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            print(f"\nERROR: command failed with return code {result.returncode}")
            sys.exit(1)


def count_sequences(fasta_path):
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Data paths
    parser.add_argument('--previous-data', required=True, help='Previous pipeline data directory (e.g. data/20260915)')
    parser.add_argument('--outdir', default=None, help='Output data directory (default: data/YYYYMMDD)')
    parser.add_argument('--prefix', default='Ituri2026.DRC_trimmed', help='Output prefix for pipeline files')

    # Exclude lists (canonical, persist across datasets)
    parser.add_argument('--exclude-identity', default='data/exclude_identity.txt',
                        help='Identity exclusions (non-consortium duplicates)')
    parser.add_argument('--exclude-quality', default='data/exclude_quality.txt',
                        help='Quality exclusions (garbage sequences)')
    parser.add_argument('--exclude-pending', default='data/exclude_pending.txt',
                        help='Pending review exclusions (temporarily excluded until confirmed or cleared)')

    # Pipeline options
    parser.add_argument('--adar-mode', default='mask', choices=['mask', 'exclude_tc', 'exclude_all'])
    parser.add_argument('--outlier-threshold', type=float, default=3.0, help='Auto-remove outliers at this z-score (default: 3.0)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--aliases-file', default='data/aliases.csv')

    # XML generation
    parser.add_argument('--template-sg', default=None, help='SkyGrid template XML')
    parser.add_argument('--template-egc', default=None, help='EGC template XML')
    parser.add_argument('--anchor-date', default='2025-11-08', help='SkyGrid anchor date')
    parser.add_argument('--chain-length', type=int, default=100_000_000)
    parser.add_argument('--log-every', type=int, default=10_000)

    # Control
    parser.add_argument('--skip-download', action='store_true', help='Skip download (use existing files in --outdir)')
    parser.add_argument('--skip-pipeline', action='store_true', help='Skip pipeline (use existing pipeline outputs)')
    parser.add_argument('--skip-beast-xmls', action='store_true', help='Skip XML generation')
    parser.add_argument('--full-check', action='store_true', help='Run health checks on all genomes, not just new ones')
    parser.add_argument('--organism', default='ebola-bdbv')

    # Notification
    parser.add_argument('--slack-webhook', default=None, help='Slack webhook URL for notifications')

    args = parser.parse_args()

    today = datetime.now().strftime('%Y%m%d')
    if args.outdir is None:
        args.outdir = f'data/{today}'
    outdir = Path(args.outdir)

    # ------------------------------------------------------------------
    # 1. Download
    # ------------------------------------------------------------------
    if not args.skip_download:
        download_script = SCRIPTS_DIR / 'download_pathoplexus.py'
        stdout = run(
            f'{PYTHON} {download_script} -o {outdir} --organism {args.organism}',
            'Downloading from Pathoplexus',
            capture=True
        )
        print(stdout)
        paths = json.loads(stdout.strip().split('\n')[-1])
        fasta_path = paths['fasta']
        metadata_path = paths['metadata']
    else:
        fasta_files = sorted(outdir.glob(f'{args.organism}_aligned-nuc_*.fasta'))
        metadata_files = sorted(outdir.glob(f'{args.organism}_metadata_*.tsv'))
        if not fasta_files or not metadata_files:
            print(f"ERROR: no FASTA/metadata files found in {outdir}")
            sys.exit(1)
        fasta_path = str(fasta_files[-1])
        metadata_path = str(metadata_files[-1])
        print(f"Using existing files:\n  FASTA: {fasta_path}\n  Metadata: {metadata_path}")

    # ------------------------------------------------------------------
    # 2. Quick count check
    # ------------------------------------------------------------------
    n_current = count_sequences(fasta_path)
    prev_fasta = sorted(Path(args.previous_data).glob('*.fasta'))
    if prev_fasta:
        n_previous = count_sequences(str(prev_fasta[-1]))
    else:
        prev_aln = sorted(Path(args.previous_data).glob(f'{args.prefix}*.final.aln'))
        n_previous = count_sequences(str(prev_aln[-1])) if prev_aln else 0

    print(f"\nSequence count: {n_current} (previous: {n_previous})")
    if n_current == n_previous and not args.skip_download:
        print("No new genomes. Exiting.")
        return
    if n_current < n_previous:
        print(f"WARNING: fewer sequences than previous ({n_current} < {n_previous})")

    # ------------------------------------------------------------------
    # 3. Merge exclude lists
    # ------------------------------------------------------------------
    merged_exclude = outdir / 'exclude_merged.txt'
    n_excluded = merge_exclude_files(
        merged_exclude,
        args.exclude_identity, args.exclude_quality, args.exclude_pending
    )
    print(f"Merged exclude list: {n_excluded} entries → {merged_exclude}")
    for label, path in [('identity', args.exclude_identity),
                        ('quality', args.exclude_quality),
                        ('pending', args.exclude_pending)]:
        n = len(read_exclude_entries(path))
        print(f"  {label}: {n} ({path})")

    # ------------------------------------------------------------------
    # 4. Pipeline
    # ------------------------------------------------------------------
    if not args.skip_pipeline:
        pipeline_script = SCRIPTS_DIR / 'run_pipeline.py'
        cmd = (
            f'{PYTHON} {pipeline_script}'
            f' -f {fasta_path}'
            f' -m {metadata_path}'
            f' -p {args.prefix}'
            f' -d {outdir}'
            f' -a {args.aliases_file}'
            f' --adar-mode {args.adar_mode}'
            f' --seed {args.seed}'
            f' --outlier-removal-threshold {args.outlier_threshold}'
            f' --exclude-file {merged_exclude}'
        )
        run(cmd, 'Running pipeline (clean → final)')

    # ------------------------------------------------------------------
    # 5. Health checks
    # ------------------------------------------------------------------
    health_script = SCRIPTS_DIR / 'health_checks.py'
    health_report = outdir / 'health_check_report.txt'
    full_flag = ' --full' if args.full_check else ''
    run(
        f'{PYTHON} {health_script} --current {outdir} --previous {args.previous_data} --prefix {args.prefix}{full_flag} -o {health_report}',
        'Running health checks' + (' (full)' if args.full_check else '')
    )

    # Carry forward pending exclusions: if health checks flagged duplicate_sample
    # on new genomes, append them to exclude_pending.txt
    pending_path = Path(args.exclude_pending)
    existing_pending = read_exclude_entries(pending_path)
    # Import health_checks to get structured results
    sys.path.insert(0, str(SCRIPTS_DIR))
    from health_checks import run_all_checks
    results = run_all_checks(str(outdir), args.previous_data, args.prefix, full=args.full_check)
    new_pending = []
    for entry in results.get('flagged_sequences', []):
        flag_names = {f['flag'] for f in entry['flags']}
        if 'duplicate_sample' in flag_names and entry['sample_id'] not in existing_pending:
            new_pending.append(entry['sample_id'])
    if new_pending:
        with open(pending_path, 'a') as f:
            for sid in new_pending:
                f.write(f'{sid}\n')
        print(f"Added {len(new_pending)} sequences to {pending_path}: {', '.join(new_pending)}")

    # ------------------------------------------------------------------
    # 6. Generate BEAST XMLs
    # ------------------------------------------------------------------
    if not args.skip_beast_xmls and any([args.template_sg, args.template_egc]):
        adar_tag = 'ADAR_masked' if args.adar_mode == 'mask' else 'ADAR_excluded'
        final_aln = outdir / f'{args.prefix}.{adar_tag}.final.aln'
        if not final_aln.exists():
            print(f"WARNING: {final_aln} not found, skipping XML generation")
        else:
            n_final = count_sequences(str(final_aln))
            beast_prefix = f'{args.prefix}.{adar_tag}.n{n_final}'

            xml_script = SCRIPTS_DIR / 'generate_beast_xml.py'
            cmd = (
                f'{PYTHON} {xml_script}'
                f' -a {final_aln}'
                f' -o {outdir}'
                f' -p {beast_prefix}'
                f' --chain-length {args.chain_length}'
                f' --log-every {args.log_every}'
            )
            if args.template_sg:
                cmd += f' --template-sg {args.template_sg}'
                cmd += f' --anchor-date {args.anchor_date}'
            if args.template_egc:
                cmd += f' --template-egc {args.template_egc}'
            run(cmd, f'Generating BEAST XMLs (n={n_final})')

            print(f"\nBEAST XMLs in: {outdir}")
            for xml in sorted(outdir.glob('*.xml')):
                print(f"  {xml}")

    # ------------------------------------------------------------------
    # 7. Summary
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  Update complete")
    print(f"{'='*60}")
    print(f"Data directory:    {outdir}")
    print(f"Health report:     {health_report}")
    if not args.skip_beast_xmls and any([args.template_sg, args.template_egc]):
        print(f"BEAST XMLs:        {outdir}/*.xml")

    # Slack notification
    if args.slack_webhook:
        try:
            from urllib.request import urlopen, Request
            report_text = health_report.read_text() if health_report.exists() else 'No report available'
            summary_lines = [l for l in report_text.split('\n') if l.strip() and not l.startswith('=')]
            payload = json.dumps({
                'text': f'Pipeline update complete: {outdir}\n```\n' + '\n'.join(summary_lines[:20]) + '\n```'
            })
            req = Request(args.slack_webhook, data=payload.encode(), headers={'Content-Type': 'application/json'})
            urlopen(req, timeout=10)
            print("Slack notification sent.")
        except Exception as e:
            print(f"Slack notification failed: {e}", file=sys.stderr)


if __name__ == '__main__':
    main()
