#!/usr/bin/env python3
"""Single entry point: download → quality report → pipeline → health checks.

Optional flags:
  --previous-data   Enables count comparison, health checks, and run diff
  --linelist        Enables duplicate resolution and linelist mismatch checks

Usage:
    python bin/run.py --outdir data/20261029 --prefix BDBV2026.outbreak --scope outbreak
    python bin/run.py --previous-data data/20261027 --outdir data/20261029 --linelist data/linelist.xlsx
"""

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent
PYTHON = sys.executable


class RunLog:
    """Append-only run log written incrementally to disk."""

    def __init__(self, path, args):
        self.path = Path(path)
        self.lines = []
        self._write(f'PIPELINE RUN LOG')
        self._write(f'Started: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
        self._write(f'Output:  {args.outdir}')
        self._write(f'Scope:   {args.scope}')
        self._write(f'Prefix:  {args.prefix}')
        self._write(f'Previous: {args.previous_data or "(none)"}')
        if args.trim_start or args.trim_end:
            t_start = args.trim_start or 1
            t_end = args.trim_end or 'end'
            self._write(f'Trim:    {t_start}–{t_end}')
        if hasattr(args, 'model') and args.model:
            self._write(f'Model:   {args.model}')
        if args.min_date:
            self._write(f'Min date: {args.min_date}')
        if hasattr(args, 'mask_sites') and args.mask_sites:
            self._write(f'Mask sites: {args.mask_sites}')
        if hasattr(args, 'mask_sites_file') and args.mask_sites_file:
            self._write(f'Mask sites file: {args.mask_sites_file}')
        self._write('')

    def step(self, title):
        self._write(f'── {title} ── [{datetime.now().strftime("%H:%M:%S")}]')

    def log(self, msg):
        self._write(f'  {msg}')

    def blank(self):
        self._write('')

    def finish(self):
        self._write('')
        self._write(f'Finished: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

    def _write(self, line):
        self.lines.append(line)
        with open(self.path, 'a') as f:
            f.write(line + '\n')


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
        f.write('# Merged exclude list (manual + auto)\n')
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


def run_streaming(cmd, rlog, description=None):
    """Run a command, stream stdout to terminal, and log step banners + key lines."""
    import re as _re
    if description:
        print(f"\n{'='*60}")
        print(f"  {description}")
        print(f"{'='*60}")
    print(f"$ {cmd}\n")
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    prev_line = ''
    for line in proc.stdout:
        print(line, end='')
        stripped = line.strip()
        # Step banner: the description line between === bars
        if prev_line.startswith('====') and stripped and not stripped.startswith('====') and not stripped.startswith('$'):
            rlog.log(f'{stripped} [{datetime.now().strftime("%H:%M:%S")}]')
        # Key stats lines
        if stripped.startswith('Cleaned:') or stripped.startswith('Filtered:'):
            rlog.log(f'  {stripped}')
        elif stripped.startswith('After ADAR:') or stripped.startswith('Auto-removal:'):
            rlog.log(f'  {stripped}')
        elif stripped.startswith('Best-fit model:') or stripped.startswith('Model:'):
            rlog.log(f'  {stripped}')
        elif stripped.startswith('Final:') and 'sequences' in stripped:
            rlog.log(f'  {stripped}')
        elif stripped.startswith('Pre-filter:'):
            rlog.log(f'  {stripped}')
        elif _re.match(r'Outliers detected:', stripped) or _re.match(r'Auto-removal:', stripped):
            rlog.log(f'  {stripped}')
        prev_line = stripped
    proc.wait()
    if proc.returncode != 0:
        print(f"\nERROR: command failed with return code {proc.returncode}")
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
    parser.add_argument('--previous-data', default=None, help='Previous pipeline data directory (enables comparison, health checks)')
    parser.add_argument('--outdir', default=None, help='Output data directory (default: data/YYYYMMDD)')
    parser.add_argument('--prefix', default='Ituri2026.DRC_trimmed', help='Output prefix for pipeline files')

    # Exclude list
    parser.add_argument('--exclude', default='data/exclude.txt',
                        help='Manual exclusion list (curated by hand)')

    # Metadata quality report
    parser.add_argument('--linelist', default=None,
                        help='INRB linelist Excel for duplicate resolution (enables metadata quality report)')
    parser.add_argument('--hz-lookup', default='data/health_area_zone_province_slim.csv',
                        help='Health area/zone/province CSV (for province name detection)')

    # Pipeline options
    parser.add_argument('--adar-mode', default='mask', choices=['mask', 'exclude_tc', 'exclude_all'])
    parser.add_argument('--trim-start', type=int, default=None, help='First alignment position to keep (1-based, inclusive)')
    parser.add_argument('--trim-end', type=int, default=None, help='Last alignment position to keep (1-based, inclusive)')
    parser.add_argument('--mask-sites', type=str, default=None,
                        help='Comma-separated 1-based positions to mask to N (e.g. 5246,5249,5251)')
    parser.add_argument('--mask-sites-file', type=str, default=None,
                        help='File with 1-based positions to mask, one per line (# comments allowed)')
    parser.add_argument('--outlier-threshold', type=float, default=3.0, help='Auto-remove outliers at this z-score (default: 3.0)')
    parser.add_argument('--model', default=None, help='Substitution model for IQ-TREE (e.g. GTR+F+R3). Default: ModelFinder')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--aliases-file', default='data/aliases.csv')

    # Control
    parser.add_argument('--skip-download', action='store_true', help='Skip download (use existing files in --outdir)')
    parser.add_argument('--skip-pipeline', action='store_true', help='Skip pipeline (use existing pipeline outputs)')
    parser.add_argument('--full-check', action='store_true', help='Run health checks on all genomes, not just new ones')
    parser.add_argument('--organism', default='ebola-bdbv')
    parser.add_argument('--scope', default='drc',
                        help='Download scope: drc, outbreak, or comma-separated countries (default: drc)')
    parser.add_argument('--min-date', default='2025-01-01',
                        help='Earliest collection date to download (default: 2025-01-01)')

    # Notification
    parser.add_argument('--slack-webhook', default=None, help='Slack webhook URL for notifications')

    args = parser.parse_args()

    today = datetime.now().strftime('%Y%m%d')
    if args.outdir is None:
        args.outdir = f'data/{today}'
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rlog = RunLog(outdir / 'run_log.txt', args)

    # ------------------------------------------------------------------
    # 1. Download
    # ------------------------------------------------------------------
    rlog.step('Download')
    if not args.skip_download:
        download_script = SCRIPTS_DIR / 'download_pathoplexus.py'
        stdout = run(
            f'{PYTHON} {download_script} -o {outdir} --organism {args.organism}'
            f' --scope {args.scope} --min-date {args.min_date}',
            'Downloading from Pathoplexus',
            capture=True
        )
        print(stdout)
        paths = json.loads(stdout.strip().split('\n')[-1])
        fasta_path = paths['fasta']
        metadata_path = paths['metadata']
        rlog.log(f'Downloaded from Pathoplexus (scope={args.scope}, min-date={args.min_date})')
    else:
        fasta_files = sorted(outdir.glob(f'{args.organism}_aligned-nuc_*.fasta'))
        metadata_files = sorted(outdir.glob(f'{args.organism}_metadata_*.tsv'))
        if not fasta_files or not metadata_files:
            print(f"ERROR: no FASTA/metadata files found in {outdir}")
            sys.exit(1)
        fasta_path = str(fasta_files[-1])
        metadata_path = str(metadata_files[-1])
        print(f"Using existing files:\n  FASTA: {fasta_path}\n  Metadata: {metadata_path}")
        rlog.log('Skipped download — using existing files')

    # ------------------------------------------------------------------
    # 2. Quick count check
    # ------------------------------------------------------------------
    n_current = count_sequences(fasta_path)
    if args.previous_data:
        prev_fasta = sorted(Path(args.previous_data).glob('*.fasta'))
        if prev_fasta:
            n_previous = count_sequences(str(prev_fasta[-1]))
        else:
            prev_aln = sorted(Path(args.previous_data).glob(f'{args.prefix}*.final.aln'))
            n_previous = count_sequences(str(prev_aln[-1])) if prev_aln else 0

        rlog.log(f'Sequences downloaded: {n_current}  (previous run: {n_previous})')
        rlog.log(f'New sequences: {n_current - n_previous}')
        rlog.blank()

        print(f"\nSequence count: {n_current} (previous: {n_previous})")
        if n_current == n_previous and not args.skip_download:
            rlog.log('No new genomes. Exiting.')
            rlog.finish()
            print("No new genomes. Exiting.")
            return
        if n_current < n_previous:
            rlog.log(f'WARNING: fewer sequences than previous ({n_current} < {n_previous})')
            print(f"WARNING: fewer sequences than previous ({n_current} < {n_previous})")
    else:
        rlog.log(f'Sequences downloaded: {n_current}')
        rlog.blank()
        print(f"\nSequence count: {n_current}")

    # ------------------------------------------------------------------
    # 3. Metadata quality report
    # ------------------------------------------------------------------
    exclude_auto = outdir / 'exclude_auto.txt'
    report_script = SCRIPTS_DIR / 'metadata_quality_report.py'
    report_txt = outdir / 'metadata_quality_report.txt'
    non_consortium_csv = outdir / 'non_consortium_sequences.csv'

    report_cmd = (
        f'{PYTHON} {report_script}'
        f' -m {metadata_path}'
        f' -f {fasta_path}'
        f' -a {args.aliases_file}'
        f' -o {outdir / "metadata_quality_report.csv"}'
        f' --exclude-out {exclude_auto}'
        f' --report-txt {report_txt}'
        f' --non-consortium-out {non_consortium_csv}'
    )
    if args.linelist:
        report_cmd += f' -l {args.linelist}'
    if Path(args.hz_lookup).exists():
        report_cmd += f' -s {args.hz_lookup}'

    rlog.step('Metadata quality report')
    run(report_cmd, 'Metadata quality report')
    print(f"Report saved to {report_txt}")

    # Log report summary
    report_csv = outdir / 'metadata_quality_report.csv'
    if report_csv.exists():
        import pandas as pd
        rdf = pd.read_csv(report_csv)
        issues = rdf['issue'].value_counts()
        n_resolved = issues.get('duplicate_resolved', 0) // 2
        n_unresolved = issues.get('duplicate_unresolved', 0) // 2
        n_input = rdf[~rdf['issue'].isin(['duplicate_resolved'])]['accessionVersion'].nunique()
        rlog.log(f'{n_resolved} resolved duplicates, {n_unresolved} unresolved')
        rlog.log(f'{n_input} accessions needing input')
    if non_consortium_csv.exists():
        import pandas as pd
        nc = pd.read_csv(non_consortium_csv)
        rlog.log(f'{len(nc)} non-consortium sequences')
    rlog.blank()

    # ------------------------------------------------------------------
    # 4. Merge exclude lists
    # ------------------------------------------------------------------
    rlog.step('Merge exclude lists')
    merged_exclude = outdir / 'exclude_merged.txt'
    exclude_sources = [args.exclude]
    if exclude_auto.exists():
        exclude_sources.append(str(exclude_auto))
    n_excluded = merge_exclude_files(merged_exclude, *exclude_sources)
    print(f"Merged exclude list: {n_excluded} entries → {merged_exclude}")
    n_manual = len(read_exclude_entries(args.exclude))
    rlog.log(f'manual: {n_manual} ({Path(args.exclude).name})')
    print(f"  manual: {n_manual} ({args.exclude})")
    if exclude_auto.exists():
        n_auto_resolved = 0
        n_auto_temp = 0
        n_auto_quality = 0
        with open(exclude_auto) as ef:
            for line in ef:
                if line.startswith('# Resolved'):
                    section = 'resolved'
                elif line.startswith('# Unresolved'):
                    section = 'temp'
                elif line.startswith('# Sequence'):
                    section = 'quality'
                elif line.strip() and not line.startswith('#'):
                    if section == 'resolved':
                        n_auto_resolved += 1
                    elif section == 'temp':
                        n_auto_temp += 1
                    elif section == 'quality':
                        n_auto_quality += 1
        parts = [f'{n_auto_resolved} resolved duplicates']
        if n_auto_temp:
            parts.append(f'{n_auto_temp} unresolved {"duplicate" if n_auto_temp == 1 else "duplicates"} ({n_auto_temp} {"sample" if n_auto_temp == 1 else "samples"} included based on completeness)')
        if n_auto_quality:
            parts.append(f'{n_auto_quality} quality failures')
        auto_summary = ', '.join(parts)
        rlog.log(f'auto-exclude: {auto_summary}')
        print(f"  auto-exclude: {auto_summary}")
    rlog.log(f'Total excluded: {n_excluded}')
    rlog.blank()

    # ------------------------------------------------------------------
    # 5. Pipeline
    # ------------------------------------------------------------------
    rlog.step('Pipeline')
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
        if args.trim_start is not None:
            cmd += f' --trim-start {args.trim_start}'
        if args.trim_end is not None:
            cmd += f' --trim-end {args.trim_end}'
        if args.model is not None:
            cmd += f' --model {args.model}'
        if args.mask_sites is not None:
            cmd += f' --mask-sites {args.mask_sites}'
        if args.mask_sites_file is not None:
            cmd += f' --mask-sites-file {args.mask_sites_file}'
        run_streaming(cmd, rlog, 'Running pipeline (clean → final)')

        adar_tag = 'ADAR_masked' if args.adar_mode == 'mask' else 'ADAR_excluded'

        # Log ADAR summary from report file
        adar_report = outdir / f'{args.prefix}.{adar_tag}.adar_report.txt'
        if adar_report.exists():
            import re
            n_branches = 0
            n_sites = 0
            for line in adar_report.read_text().splitlines():
                if line.startswith('Branches with editing signatures:'):
                    n_branches = int(line.split(':')[1].strip())
                elif line.startswith('Raccoon entries:'):
                    m = re.search(r'(\d+)', line.split(':')[1])
                    if m:
                        n_sites = int(m.group(1))
            rlog.log(f'ADAR: {n_branches} branches, {n_sites} sites masked (see {adar_report.name})')

        # Log outliers (read from outliers.txt which has the threshold-filtered set)
        outlier_file = outdir / 'outliers.txt'
        detected_file = outdir / f'{args.prefix}.{adar_tag}.outliers.txt'
        if outlier_file.exists():
            outliers = [l.strip() for l in outlier_file.read_text().splitlines()
                        if l.strip() and not l.startswith('#')]
            n_detected = 0
            if detected_file.exists():
                n_detected = sum(1 for l in detected_file.read_text().splitlines()
                                 if l.strip() and not l.startswith('#'))
            rlog.log(f'Outliers removed: {len(outliers)} (|z| >= {args.outlier_threshold}, from {n_detected} detected at |z| >= 2)')
            for o in outliers:
                rlog.log(f'  {o}')

        # Log final count
        final_aln = outdir / f'{args.prefix}.{adar_tag}.final.aln'
        if final_aln.exists():
            n_final = count_sequences(str(final_aln))
            rlog.log(f'Final alignment: {n_final} sequences')

        # Log RTT summary
        rtt_summary = outdir / f'{args.prefix}.{adar_tag}.final.rtt_summary.txt'
        if rtt_summary.exists():
            for line in rtt_summary.read_text().splitlines():
                if line.strip():
                    rlog.log(f'RTT: {line.strip()}')
    else:
        rlog.log('Skipped pipeline')
    rlog.blank()

    # ------------------------------------------------------------------
    # 6. Health checks
    # ------------------------------------------------------------------
    rlog.step('Health checks')
    health_report = outdir / 'health_check_report.txt'
    if args.previous_data:
        health_script = SCRIPTS_DIR / 'health_checks.py'
        full_flag = ' --full' if args.full_check else ''
        run(
            f'{PYTHON} {health_script} --current {outdir} --previous {args.previous_data} --prefix {args.prefix}{full_flag} -o {health_report}',
            'Running health checks' + (' (full)' if args.full_check else '')
        )
        if health_report.exists():
            for line in health_report.read_text().splitlines():
                if line.strip() and not line.startswith('=') and not line.startswith('-'):
                    rlog.log(line.strip())
    else:
        rlog.log('Skipped (no previous run)')
    rlog.blank()

    # ------------------------------------------------------------------
    # 7. Run comparison
    # ------------------------------------------------------------------
    rlog.step('Run comparison')
    comparison_report = outdir / 'run_comparison.txt'
    if args.previous_data:
        compare_script = SCRIPTS_DIR / 'compare_runs.py'
        if not args.skip_pipeline:
            adar_tag = 'ADAR_masked' if args.adar_mode == 'mask' else 'ADAR_excluded'
        run(
            f'{PYTHON} {compare_script}'
            f' --current {outdir}'
            f' --previous {args.previous_data}'
            f' --prefix {args.prefix}'
            f' --adar-mode {args.adar_mode}'
            f' -o {comparison_report}',
            'Comparing to previous run'
        )
        if comparison_report.exists():
            for line in comparison_report.read_text().splitlines():
                if line.strip() and not line.startswith('=') and not line.startswith('-'):
                    rlog.log(line.strip())
    else:
        rlog.log('Skipped (no previous run)')
    rlog.blank()

    # ------------------------------------------------------------------
    # 8. Summary
    # ------------------------------------------------------------------
    rlog.step('Output files')
    for f in sorted(outdir.iterdir()):
        if f.suffix in ('.txt', '.csv', '.tsv', '.aln', '.fasta', '.treefile', '.pdf', '.png'):
            rlog.log(f.name)
    rlog.blank()
    rlog.finish()

    print(f"\n{'='*60}")
    print(f"  Run complete")
    print(f"{'='*60}")
    print(f"Data directory:    {outdir}")
    print(f"Run log:           {outdir / 'run_log.txt'}")
    if args.previous_data:
        print(f"Health report:     {health_report}")
        print(f"Run comparison:    {comparison_report}")
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
