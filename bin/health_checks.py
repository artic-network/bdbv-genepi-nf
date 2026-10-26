#!/usr/bin/env python3
"""Modular health checks for pipeline outputs.

Compares current and previous pipeline runs. Reports per-sequence flags
(new genomes only) and dataset-level checks.

Usage:
    python scripts/health_checks.py --current data/20261012 --previous data/20260915
"""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_file(directory, pattern):
    matches = sorted(Path(directory).rglob(pattern))
    return matches[-1] if matches else None


def parse_fasta_headers(fasta_path):
    headers = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:].strip())
    return headers


def parse_fasta_sequences(fasta_path):
    seqs = {}
    current = None
    parts = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if current is not None:
                    seqs[current] = ''.join(parts)
                current = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
    if current is not None:
        seqs[current] = ''.join(parts)
    return seqs


def parse_rtt_summary(path):
    vals = {}
    with open(path) as f:
        for line in f:
            if line.startswith('Clock rate:'):
                vals['rate'] = float(line.split(':')[1].strip().split()[0])
            elif line.startswith('R-squared:'):
                vals['r_squared'] = float(line.split(':')[1].strip())
            elif line.startswith('Tips:'):
                vals['tips'] = int(line.split(':')[1].strip())
            elif line.startswith('Outliers'):
                m = re.search(r':\s*(\d+)', line)
                if m:
                    vals['outlier_count'] = int(m.group(1))
    return vals


def parse_root_to_tip(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows


def parse_outliers(path):
    outliers = {}
    if not path.exists():
        return outliers
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            m = re.search(r'z=([+-]?[\d.]+)', line)
            if m:
                name = line.split('\t')[0].strip()
                outliers[name] = float(m.group(1))
    return outliers


def parse_adar_report(path):
    info = {'n_branches': 0, 'branch_details': {}}
    if not path or not Path(path).exists():
        return info
    in_details = False
    with open(path) as f:
        for line in f:
            if line.startswith('Branches with editing signatures:'):
                info['n_branches'] = int(line.split(':')[1].strip())
            elif line.startswith('Per-branch details:'):
                in_details = True
            elif in_details and line.strip().startswith(('BIA-', 'MG-', 'NIZI-', 'NIAN-', '26FHV', 'NIZ-', 'Node')):
                m = re.match(r'\s+(\S+)\s+\((\d+)\s+T>C\)', line)
                if m:
                    info['branch_details'][m.group(1)] = int(m.group(2))
    return info


def accession_from_header(header):
    parts = header.split('|')
    if len(parts) > 1:
        return parts[1].split('.')[0], parts[1]
    return header.split('.')[0], header


def sample_id_from_header(header):
    return header.split('|')[0]


# ---------------------------------------------------------------------------
# Per-sequence checks (new genomes only)
# ---------------------------------------------------------------------------

SEQUENCE_CHECKS = []

def sequence_check(fn):
    SEQUENCE_CHECKS.append(fn)
    return fn


@sequence_check
def check_high_adar(header, ctx):
    adar = ctx.get('adar_report', {})
    sample_id = sample_id_from_header(header)
    count = adar.get('branch_details', {}).get(sample_id, 0)
    if count >= 5:
        return {'flag': 'high_adar', 'value': count, 'detail': f'{count} T>C edits on branch'}
    return None


@sequence_check
def check_high_Ns(header, ctx):
    seq = ctx.get('sequences', {}).get(header, '')
    if not seq:
        return None
    aln_len = len(seq)
    n_count = seq.upper().count('N')
    frac = n_count / aln_len if aln_len > 0 else 0
    if frac > 0.10:
        return {'flag': 'high_Ns', 'value': f'{frac:.1%}', 'detail': f'{n_count}/{aln_len} positions are N'}
    return None


@sequence_check
def check_low_coverage(header, ctx):
    seq = ctx.get('sequences', {}).get(header, '')
    if not seq:
        return None
    aln_len = len(seq)
    called = sum(1 for c in seq.upper() if c in 'ACGT')
    cov = called / aln_len if aln_len > 0 else 0
    if cov < 0.80:
        return {'flag': 'low_coverage', 'value': f'{cov:.1%}', 'detail': f'{called}/{aln_len} positions called'}
    return None


@sequence_check
def check_near_outlier(header, ctx):
    rtt = ctx.get('rtt_by_name', {})
    row = rtt.get(header)
    if row is None:
        return None
    z = abs(float(row.get('residual z-score', 0)))
    if 2.5 < z < 3.0:
        return {'flag': 'near_outlier', 'value': f'z={z:.2f}', 'detail': 'Close to removal threshold (|z| between 2.5 and 3.0)'}
    return None


@sequence_check
def check_many_gaps(header, ctx):
    seq = ctx.get('sequences', {}).get(header, '')
    if not seq:
        return None
    gap_count = seq.count('-')
    if gap_count > 500:
        return {'flag': 'many_gaps', 'value': gap_count, 'detail': f'{gap_count} gap positions'}
    return None


@sequence_check
def check_long_branch(header, ctx):
    rtt = ctx.get('rtt_by_name', {})
    row = rtt.get(header)
    if row is None:
        return None
    div = float(row.get('divergence', 0))
    median_div = ctx.get('median_divergence', 0)
    if median_div > 0 and div > 3 * median_div:
        return {'flag': 'long_branch', 'value': f'{div:.6f}', 'detail': f'Divergence {div:.6f} > 3× median {median_div:.6f}'}
    return None


@sequence_check
def check_duplicate_sample(header, ctx):
    sample_id = sample_id_from_header(header)
    dupes = ctx.get('sample_id_counts', {})
    if dupes.get(sample_id, 0) > 1:
        return {'flag': 'duplicate_sample', 'value': dupes[sample_id], 'detail': f'Sample ID {sample_id} appears {dupes[sample_id]} times'}
    return None


@sequence_check
def check_phylo_conflict(header, ctx):
    phylo_flags = ctx.get('phylo_flag_seqs', set())
    sample_id = sample_id_from_header(header)
    if sample_id in phylo_flags:
        return {'flag': 'phylo_conflict', 'value': True, 'detail': 'Flagged by raccoon phylo_flags (long_branch or high homoplasy)'}
    return None


# ---------------------------------------------------------------------------
# Dataset-level checks
# ---------------------------------------------------------------------------

DATASET_CHECKS = []

def dataset_check(fn):
    DATASET_CHECKS.append(fn)
    return fn


@dataset_check
def check_outlier_count(current, previous):
    cur = current.get('outlier_count', 0)
    prev = previous.get('outlier_count', 0)
    n_new = current.get('n_new_genomes', 0)

    if cur > 5:
        n_new_outliers = cur - prev if prev else cur
        if n_new > 0 and n_new_outliers > max(2, n_new * 0.05):
            return ('alert', f'{cur} total outliers ({n_new_outliers} new among {n_new} new genomes — unexpectedly many)')
        return ('warn', f'{cur} total outliers (was {prev})')
    return ('ok', f'{cur} outliers')


@dataset_check
def check_rate_shift(current, previous):
    cur_rate = current.get('rate')
    prev_rate = previous.get('rate')
    if cur_rate is None:
        return ('warn', 'No clock rate available')
    if prev_rate is None:
        return ('ok', f'Clock rate: {cur_rate:.2e} (no previous to compare)')
    change = abs(cur_rate - prev_rate) / prev_rate
    if change > 0.30:
        return ('alert', f'Clock rate changed {change:.0%}: {prev_rate:.2e} → {cur_rate:.2e}')
    return ('ok', f'Clock rate: {cur_rate:.2e} (Δ{change:.0%} from {prev_rate:.2e})')


@dataset_check
def check_r_squared(current, previous):
    cur_r2 = current.get('r_squared')
    prev_r2 = previous.get('r_squared')
    if cur_r2 is None:
        return ('warn', 'No R² available')
    if cur_r2 < 0.3:
        return ('warn', f'R² = {cur_r2:.3f} (below 0.3)')
    if prev_r2 is not None and (prev_r2 - cur_r2) > 0.1:
        return ('warn', f'R² dropped {prev_r2:.3f} → {cur_r2:.3f}')
    return ('ok', f'R² = {cur_r2:.3f}')


@dataset_check
def check_adar_count(current, previous):
    cur_n = current.get('adar_branches', 0)
    prev_n = previous.get('adar_branches', 0)
    if prev_n > 0 and cur_n > prev_n * 1.5:
        return ('warn', f'ADAR branches increased {prev_n} → {cur_n} (+{cur_n - prev_n})')
    return ('ok', f'{cur_n} ADAR branches')


@dataset_check
def check_genome_count(current, previous):
    cur_n = current.get('tips', 0)
    prev_n = previous.get('tips', 0)
    n_new = cur_n - prev_n
    if cur_n < prev_n:
        return ('alert', f'Fewer genomes: {prev_n} → {cur_n} (lost {prev_n - cur_n})')
    if n_new > 200:
        return ('warn', f'{n_new} new genomes in one batch ({prev_n} → {cur_n})')
    return ('ok', f'{cur_n} genomes (+{n_new} new)')


@dataset_check
def check_version_changes(current, previous):
    removed = current.get('removed_accessions', [])
    updated = current.get('updated_accessions', [])
    parts = []
    status = 'ok'
    if removed:
        status = 'alert'
        parts.append(f'{len(removed)} retracted: {", ".join(removed[:5])}')
    if updated:
        status = 'warn' if status == 'ok' else status
        parts.append(f'{len(updated)} version-updated')
    if not parts:
        return ('ok', 'No retractions or version changes')
    return (status, '; '.join(parts))


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------

def load_context(data_dir, prefix_pattern='Ituri2026.DRC_trimmed'):
    data_dir = Path(data_dir)
    ctx = {}

    rtt_file = find_file(data_dir, f'{prefix_pattern}.ADAR_masked.bestFitRoot.rootToTip.csv')
    if rtt_file:
        rtt_rows = parse_root_to_tip(rtt_file)
        ctx['rtt_by_name'] = {row['name']: row for row in rtt_rows}
        divs = [float(row['divergence']) for row in rtt_rows if row.get('divergence')]
        if divs:
            divs.sort()
            ctx['median_divergence'] = divs[len(divs) // 2]

    rtt_summary = find_file(data_dir, f'{prefix_pattern}.ADAR_masked.rtt_summary.txt')
    if rtt_summary:
        ctx['rtt_summary'] = parse_rtt_summary(rtt_summary)

    final_rtt_summary = find_file(data_dir, f'{prefix_pattern}.ADAR_masked.final.rtt_summary.txt')
    if final_rtt_summary:
        ctx['final_rtt_summary'] = parse_rtt_summary(final_rtt_summary)

    adar_file = find_file(data_dir, f'{prefix_pattern}.ADAR_masked.adar_report.txt')
    if adar_file:
        ctx['adar_report'] = parse_adar_report(adar_file)

    phylo_flags_file = find_file(data_dir, '*.phylo_flags.csv')
    if phylo_flags_file:
        phylo_seqs = set()
        with open(phylo_flags_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get('mask_boolean', '').lower() == 'true' and 'long_branch' in row.get('mutation_type', ''):
                    for entry in row.get('present_in', '').split(';'):
                        sample = entry.split('|')[0].strip()
                        if sample and not sample.startswith('Node'):
                            phylo_seqs.add(sample)
        ctx['phylo_flag_seqs'] = phylo_seqs

    outlier_file = find_file(data_dir, f'{prefix_pattern}.ADAR_masked.outliers.txt')
    if outlier_file:
        ctx['outliers'] = parse_outliers(outlier_file)

    fasta_file = find_file(data_dir, f'{prefix_pattern}.ADAR_masked.aln')
    if not fasta_file:
        fasta_file = find_file(data_dir, f'{prefix_pattern}.aln')
    if fasta_file:
        ctx['sequences'] = parse_fasta_sequences(fasta_file)
        headers = list(ctx['sequences'].keys())
        sample_counts = defaultdict(int)
        for h in headers:
            sample_counts[sample_id_from_header(h)] += 1
        ctx['sample_id_counts'] = dict(sample_counts)

    metadata_file = find_file(data_dir, '*_metadata_*.tsv')
    if metadata_file:
        accessions = {}
        with open(metadata_file) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                av = row.get('accessionVersion', '')
                if av:
                    acc = av.rsplit('.', 1)[0]
                    ver = av.rsplit('.', 1)[1] if '.' in av else '1'
                    accessions[acc] = ver
        ctx['accessions'] = accessions

    return ctx


def load_watchlist(data_dir):
    watchlist_file = Path(data_dir).parent / 'watchlist.txt'
    if not watchlist_file.exists():
        return {}
    entries = {}
    for line in watchlist_file.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            entries[parts[0]] = parts[1]
    return entries


def run_all_checks(current_dir, previous_dir, prefix_pattern='Ituri2026.DRC_trimmed', full=False):
    cur_ctx = load_context(current_dir, prefix_pattern)
    prev_ctx = load_context(previous_dir, prefix_pattern) if previous_dir else {}

    # Determine which genomes to flag
    cur_headers = set(cur_ctx.get('sequences', {}).keys())
    prev_headers = set(prev_ctx.get('sequences', {}).keys())
    new_headers = cur_headers - prev_headers
    check_headers = cur_headers if full else new_headers

    # Version tracking
    cur_acc = cur_ctx.get('accessions', {})
    prev_acc = prev_ctx.get('accessions', {})
    added = sorted(set(cur_acc.keys()) - set(prev_acc.keys()))
    removed = sorted(set(prev_acc.keys()) - set(cur_acc.keys()))
    updated = sorted(a for a in set(cur_acc.keys()) & set(prev_acc.keys()) if cur_acc[a] != prev_acc[a])

    # Per-sequence flags
    flagged = []
    for header in sorted(check_headers):
        flags = []
        for check_fn in SEQUENCE_CHECKS:
            result = check_fn(header, cur_ctx)
            if result:
                flags.append(result)
        if flags:
            flagged.append({'header': header, 'sample_id': sample_id_from_header(header), 'flags': flags})

    # Dataset-level checks
    cur_summary = cur_ctx.get('rtt_summary', {})
    prev_summary = prev_ctx.get('rtt_summary', {})
    cur_summary['adar_branches'] = cur_ctx.get('adar_report', {}).get('n_branches', 0)
    prev_summary['adar_branches'] = prev_ctx.get('adar_report', {}).get('n_branches', 0)
    cur_summary['outlier_count'] = len(cur_ctx.get('outliers', {}))
    prev_summary['outlier_count'] = len(prev_ctx.get('outliers', {}))
    cur_summary['n_new_genomes'] = len(new_headers)
    cur_summary['removed_accessions'] = removed
    cur_summary['updated_accessions'] = updated

    dataset_results = []
    for check_fn in DATASET_CHECKS:
        name = check_fn.__name__.replace('check_', '')
        status, message = check_fn(cur_summary, prev_summary)
        dataset_results.append((name, status, message))

    watchlist = load_watchlist(current_dir)
    watchlist_status = []
    cur_rtt = cur_ctx.get('rtt_by_name', {})
    prev_rtt = prev_ctx.get('rtt_by_name', {}) if prev_ctx else {}
    for acc, sample_id in watchlist.items():
        cur_match = [k for k in cur_rtt if acc in k]
        if not cur_match:
            watchlist_status.append((sample_id, acc, 'missing', 'Not in current dataset'))
            continue
        cur_row = cur_rtt[cur_match[0]]
        cur_z = float(cur_row.get('residual z-score', cur_row.get('z_score', 0)))
        prev_match = [k for k in prev_rtt if acc in k]
        if prev_match:
            prev_z = float(prev_rtt[prev_match[0]].get('residual z-score', prev_rtt[prev_match[0]].get('z_score', 0)))
            delta = cur_z - prev_z
            status_str = f'z={cur_z:+.2f} (was {prev_z:+.2f}, Δ{delta:+.2f})'
        else:
            status_str = f'z={cur_z:+.2f} (new to comparison)'
        if abs(cur_z) >= 3:
            watchlist_status.append((sample_id, acc, 'alert', status_str))
        elif abs(cur_z) >= 2.5:
            watchlist_status.append((sample_id, acc, 'warn', status_str))
        else:
            watchlist_status.append((sample_id, acc, 'ok', status_str))

    return {
        'flagged_sequences': flagged,
        'version_changes': {'added': added, 'removed': removed, 'updated': updated},
        'dataset_checks': dataset_results,
        'watchlist': watchlist_status,
        'n_new': len(new_headers),
        'n_checked': len(check_headers),
        'n_total': len(cur_headers),
        'full': full,
    }


def format_report(results):
    lines = []
    lines.append('=' * 60)
    lines.append('  Health Check Report')
    lines.append('=' * 60)
    if results.get('full'):
        lines.append(f'\nGenomes: {results["n_total"]} total (full check)\n')
    else:
        lines.append(f'\nGenomes: {results["n_total"]} total, {results["n_new"]} new\n')

    lines.append('--- Dataset checks ---')
    max_severity = 'ok'
    for name, status, message in results['dataset_checks']:
        icon = {'ok': ' OK ', 'warn': 'WARN', 'alert': '!!!!'}[status]
        lines.append(f'  [{icon}] {name}: {message}')
        if status == 'alert':
            max_severity = 'alert'
        elif status == 'warn' and max_severity != 'alert':
            max_severity = 'warn'

    vc = results['version_changes']
    lines.append(f'\n--- Version tracking ---')
    lines.append(f'  Added: {len(vc["added"])} new accessions')
    if vc['removed']:
        lines.append(f'  RETRACTED: {", ".join(vc["removed"])}')
    if vc['updated']:
        lines.append(f'  Updated versions: {", ".join(vc["updated"])}')

    flagged = results['flagged_sequences']
    n_checked = results.get('n_checked', results['n_new'])
    scope = 'all' if results.get('full') else 'new'
    lines.append(f'\n--- Flagged sequences ({len(flagged)}/{n_checked} {scope}) ---')
    if flagged:
        flag_counts = defaultdict(int)
        for entry in flagged:
            for f in entry['flags']:
                flag_counts[f['flag']] += 1
        lines.append(f'  Flag summary: {", ".join(f"{k}={v}" for k, v in sorted(flag_counts.items(), key=lambda x: -x[1]))}')
        lines.append('')
        for entry in flagged:
            flag_strs = [f'{f["flag"]}({f["value"]})' for f in entry['flags']]
            lines.append(f'  {entry["sample_id"]}: {", ".join(flag_strs)}')
    else:
        lines.append('  None flagged.')

    watchlist = results.get('watchlist', [])
    if watchlist:
        lines.append(f'\n--- Watchlist ({len(watchlist)}) ---')
        for sample_id, acc, status, detail in watchlist:
            icon = {'ok': ' OK ', 'warn': 'WARN', 'alert': '!!!!', 'missing': ' ?? '}[status]
            lines.append(f'  [{icon}] {sample_id} ({acc}): {detail}')
            if status == 'alert':
                max_severity = 'alert'

    lines.append(f'\nOverall status: {max_severity.upper()}')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--current', required=True, help='Current pipeline output directory')
    parser.add_argument('--previous', default=None, help='Previous pipeline output directory (for comparison)')
    parser.add_argument('--prefix', default='Ituri2026.DRC_trimmed', help='File prefix pattern')
    parser.add_argument('--full', action='store_true',
                        help='Check all genomes, not just new ones since previous')
    parser.add_argument('-o', '--output', default=None, help='Save report to file')
    args = parser.parse_args()

    results = run_all_checks(args.current, args.previous, args.prefix, full=args.full)
    report = format_report(results)
    print(report)

    if args.output:
        Path(args.output).write_text(report + '\n')
        print(f'\nReport saved to {args.output}')

    retracted = results['version_changes'].get('removed', [])
    if retracted:
        retractions_file = Path(args.current).parent / 'retractions.txt'
        existing = set()
        if retractions_file.exists():
            for line in retractions_file.read_text().splitlines():
                if line and not line.startswith('#'):
                    existing.add(line.split('\t')[0])
        new_retractions = [a for a in retracted if a not in existing]
        if new_retractions:
            today = Path(args.current).name
            with open(retractions_file, 'a') as f:
                for acc in new_retractions:
                    f.write(f'{acc}\t{today}\t\n')


if __name__ == '__main__':
    main()
