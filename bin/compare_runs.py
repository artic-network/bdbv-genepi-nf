#!/usr/bin/env python3
"""Compare two pipeline runs side-by-side.

Produces a human-readable summary of what changed between runs:
sequence counts, added/dropped accessions, exclude list changes,
clock signal, ADAR edits, metadata updates, and quality report diffs.

Usage:
    python bin/compare_runs.py \
        --current data/20261027 \
        --previous data/20261012 \
        --prefix Ituri2026.DRC_trimmed \
        -o data/20261027/run_comparison.txt
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def count_sequences(path):
    if not path or not Path(path).exists():
        return None
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                n += 1
    return n


def get_accessions(path):
    if not path or not Path(path).exists():
        return set()
    accs = set()
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].strip().split('|')
                accs.add(parts[1] if len(parts) > 1 else parts[0])
    return accs


def read_exclude(path):
    entries = set()
    if not path or not Path(path).exists():
        return entries
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                entries.add(line.split()[0])
    return entries


def parse_rtt_summary(path):
    if not path or not Path(path).exists():
        return {}
    vals = {}
    with open(path) as f:
        for line in f:
            if ':' in line:
                key, val = line.split(':', 1)
                key = key.strip()
                val = val.strip()
                if key == 'Clock rate':
                    vals['clock_rate'] = val.split()[0]
                elif key == 'R-squared':
                    vals['r_squared'] = val
                elif key == 'Root date':
                    vals['root_date'] = val
                elif key == 'Tips':
                    vals['tips'] = val
                elif key.startswith('Outliers'):
                    vals['outliers'] = val
    return vals


def parse_adar_summary(path):
    if not path or not Path(path).exists():
        return {}
    vals = {}
    with open(path) as f:
        for line in f:
            if line.startswith('Input:'):
                vals['input'] = line.split(':')[1].strip().split()[0]
            elif line.startswith('Branches with editing signatures:'):
                vals['branches'] = line.split(':')[1].strip()
            elif line.startswith('Raccoon entries:'):
                vals['raccoon_entries'] = line.split(':')[1].strip().split()[0]
    return vals


def find_file(directory, patterns):
    d = Path(directory)
    for pat in patterns:
        matches = sorted(d.glob(pat))
        if matches:
            return matches[-1]
    return None


def find_metadata(directory):
    d = Path(directory)
    for pat in ['metadata.tsv', 'ebola-bdbv_metadata_*.tsv', '*_metadata_*.tsv']:
        matches = sorted(d.glob(pat))
        if matches:
            return matches[-1]
    return None


def fmt_count(n):
    if n is None:
        return '—'
    return str(n)


def fmt_delta(prev, curr):
    if prev is None or curr is None:
        return ''
    d = curr - prev
    if d > 0:
        return f'(+{d})'
    elif d < 0:
        return f'({d})'
    return '(=)'


def compare(current_dir, previous_dir, prefix, adar_tag, output_path):
    lines = []

    def pr(s=''):
        lines.append(s)

    cur = Path(current_dir)

    if not previous_dir or not Path(previous_dir).exists():
        pr('No previous run to compare against. Skipping comparison.')
        report = '\n'.join(lines)
        print(report)
        if output_path:
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(report + '\n')
        return

    prev = Path(previous_dir)

    pr(f'Run comparison: {prev.name} → {cur.name}')
    pr('=' * 60)

    # --- 1. Sequence counts ---
    stages = [
        ('Downloaded', ['*.fasta', 'aligned.fasta']),
        ('Cleaned', [f'{prefix}.cleaned.aln']),
        ('After exclusions', [f'{prefix}.aln']),
        ('After ADAR', [f'{prefix}.{adar_tag}.aln']),
        ('Final', [f'{prefix}.{adar_tag}.final.aln']),
    ]

    pr()
    pr('1. SEQUENCE COUNTS')
    pr('-' * 40)
    pr(f'  {"Stage":<22} {"Previous":>10} {"Current":>10} {"Change":>10}')
    pr(f'  {"—"*22} {"—"*10} {"—"*10} {"—"*10}')

    for label, pats in stages:
        pf = find_file(prev, pats)
        cf = find_file(cur, pats)
        np = count_sequences(pf) if pf else None
        nc = count_sequences(cf) if cf else None
        pr(f'  {label:<22} {fmt_count(np):>10} {fmt_count(nc):>10} {fmt_delta(np, nc):>10}')

    # --- 2. Sequence-level diff ---
    prev_final = find_file(prev, [f'{prefix}.{adar_tag}.final.aln'])
    cur_final = find_file(cur, [f'{prefix}.{adar_tag}.final.aln'])

    prev_accs = get_accessions(prev_final) if prev_final else set()
    cur_accs = get_accessions(cur_final) if cur_final else set()

    added = sorted(cur_accs - prev_accs)
    dropped = sorted(prev_accs - cur_accs)

    cur_excludes = read_exclude(cur / 'exclude_merged.txt')
    cur_outliers = set()
    outlier_file = find_file(cur, ['outliers.txt', f'{prefix}.{adar_tag}.outliers.txt'])
    if outlier_file:
        with open(outlier_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('|')
                    if len(parts) > 1:
                        cur_outliers.add(parts[1])
                    else:
                        cur_outliers.add(parts[0].split('\t')[0])

    pr()
    pr('2. SEQUENCE-LEVEL CHANGES')
    pr('-' * 40)
    pr(f'  Added:   {len(added)}')
    pr(f'  Dropped: {len(dropped)}')

    if added:
        pr()
        pr('  New accessions:')
        for acc in added[:50]:
            pr(f'    + {acc}')
        if len(added) > 50:
            pr(f'    ... and {len(added) - 50} more')

    if dropped:
        pr()
        pr('  Dropped accessions:')
        for acc in dropped:
            reason = ''
            if acc in cur_excludes:
                reason = ' [excluded]'
            elif acc in cur_outliers:
                reason = ' [outlier]'
            pr(f'    - {acc}{reason}')

    # --- 3. Exclude list diff ---
    pr()
    pr('3. EXCLUDE LIST CHANGES')
    pr('-' * 40)

    prev_excl = read_exclude(prev / 'exclude_merged.txt')
    cur_excl = read_exclude(cur / 'exclude_merged.txt')

    pr(f'  Previous: {len(prev_excl)} entries')
    pr(f'  Current:  {len(cur_excl)} entries')

    new_excl = sorted(cur_excl - prev_excl)
    removed_excl = sorted(prev_excl - cur_excl)

    if new_excl:
        pr(f'  Newly excluded ({len(new_excl)}):')
        for e in new_excl:
            pr(f'    + {e}')
    if removed_excl:
        pr(f'  No longer excluded ({len(removed_excl)}):')
        for e in removed_excl:
            pr(f'    - {e}')
    if not new_excl and not removed_excl:
        pr('  No changes')

    # --- 4. Clock signal ---
    pr()
    pr('4. CLOCK SIGNAL')
    pr('-' * 40)

    prev_rtt = parse_rtt_summary(
        find_file(prev, [f'{prefix}.{adar_tag}.final.rtt_summary.txt']) or ''
    )
    cur_rtt = parse_rtt_summary(
        find_file(cur, [f'{prefix}.{adar_tag}.final.rtt_summary.txt']) or ''
    )

    if prev_rtt or cur_rtt:
        for key, label in [('clock_rate', 'Clock rate'),
                           ('r_squared', 'R²'),
                           ('root_date', 'Root date'),
                           ('tips', 'Tips'),
                           ('outliers', 'Outliers')]:
            pv = prev_rtt.get(key, '—')
            cv = cur_rtt.get(key, '—')
            pr(f'  {label:<15} {pv:>25} → {cv}')
    else:
        pr('  No RTT summary files found')

    # --- 5. ADAR summary ---
    pr()
    pr('5. ADAR EDITS')
    pr('-' * 40)

    prev_adar = parse_adar_summary(
        find_file(prev, [f'{prefix}.{adar_tag}.adar_report.txt']) or ''
    )
    cur_adar = parse_adar_summary(
        find_file(cur, [f'{prefix}.{adar_tag}.adar_report.txt']) or ''
    )

    if prev_adar or cur_adar:
        for key, label in [('branches', 'Branches with ADAR'),
                           ('raccoon_entries', 'Raccoon entries')]:
            pv = prev_adar.get(key, '—')
            cv = cur_adar.get(key, '—')
            pr(f'  {label:<25} {pv:>10} → {cv}')
    else:
        pr('  No ADAR report files found')

    # --- 6. Metadata changes ---
    pr()
    pr('6. METADATA CHANGES')
    pr('-' * 40)

    prev_meta_f = find_metadata(prev)
    cur_meta_f = find_metadata(cur)

    if prev_meta_f and cur_meta_f:
        pm = pd.read_csv(prev_meta_f, sep='\t', dtype=str)
        cm = pd.read_csv(cur_meta_f, sep='\t', dtype=str)

        common = set(pm['accessionVersion']) & set(cm['accessionVersion'])
        pm_idx = pm.set_index('accessionVersion')
        cm_idx = cm.set_index('accessionVersion')

        changes = []
        for acc in sorted(common):
            for col in ['geoLocAdmin2', 'sampleCollectionDate']:
                if col not in pm_idx.columns or col not in cm_idx.columns:
                    continue
                old = pm_idx.loc[acc, col] if acc in pm_idx.index else None
                new = cm_idx.loc[acc, col] if acc in cm_idx.index else None
                if str(old) != str(new):
                    changes.append((acc, col, old, new))

        if changes:
            pr(f'  {len(changes)} field(s) changed across {len(set(c[0] for c in changes))} accession(s):')
            for acc, col, old, new in changes[:30]:
                field = 'location' if 'Admin' in col else 'date'
                pr(f'    {acc}: {field} "{old}" → "{new}"')
            if len(changes) > 30:
                pr(f'    ... and {len(changes) - 30} more')
        else:
            pr('  No location or date changes in shared accessions')
    else:
        pr('  Metadata files not found in both directories')

    # --- 7. Quality report diff ---
    pr()
    pr('7. QUALITY REPORT')
    pr('-' * 40)

    prev_qr = find_file(prev, ['metadata_quality_report.csv'])
    cur_qr = find_file(cur, ['metadata_quality_report.csv'])

    if prev_qr and cur_qr:
        pq = pd.read_csv(prev_qr)
        cq = pd.read_csv(cur_qr)

        prev_counts = pq['issue'].value_counts().to_dict()
        cur_counts = cq['issue'].value_counts().to_dict()
        all_issues = sorted(set(prev_counts) | set(cur_counts))

        pr(f'  {"Issue type":<30} {"Previous":>10} {"Current":>10} {"Change":>10}')
        pr(f'  {"—"*30} {"—"*10} {"—"*10} {"—"*10}')
        for issue in all_issues:
            np = prev_counts.get(issue, 0)
            nc = cur_counts.get(issue, 0)
            pr(f'  {issue:<30} {np:>10} {nc:>10} {fmt_delta(np, nc):>10}')
        pr(f'  {"TOTAL":<30} {len(pq):>10} {len(cq):>10} {fmt_delta(len(pq), len(cq)):>10}')
    elif cur_qr:
        cq = pd.read_csv(cur_qr)
        pr(f'  Current run: {len(cq)} issues (no previous report to compare)')
        for issue, count in cq['issue'].value_counts().items():
            pr(f'    {issue}: {count}')
    else:
        pr('  No quality report CSV found')

    pr()
    pr('=' * 60)

    report = '\n'.join(lines)
    print(report)

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(report + '\n')
        print(f'\nSaved to {output_path}', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--current', required=True, help='Current run directory')
    parser.add_argument('--previous', default=None, help='Previous run directory (skips comparison if absent)')
    parser.add_argument('--prefix', default='Ituri2026.DRC_trimmed',
                        help='Alignment file prefix (default: Ituri2026.DRC_trimmed)')
    parser.add_argument('--adar-mode', default='mask',
                        choices=['mask', 'exclude_tc', 'exclude_all'])
    parser.add_argument('-o', '--output', default=None,
                        help='Output file (default: {current}/run_comparison.txt)')
    args = parser.parse_args()

    adar_tag = 'ADAR_masked' if args.adar_mode == 'mask' else 'ADAR_excluded'
    output = args.output or str(Path(args.current) / 'run_comparison.txt')

    compare(args.current, args.previous, args.prefix, adar_tag, output)


if __name__ == '__main__':
    main()
