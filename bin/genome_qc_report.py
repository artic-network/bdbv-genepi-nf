#!/usr/bin/env python3
import argparse
import csv
import re
import statistics
from gzip import open as gopen


def readFASTA(path):
    if path.endswith('.gz'):
        stream = gopen(path)
    else:
        stream = open(path)
    seqs = {}
    name = None
    seq = ''
    for line in stream:
        if isinstance(line, bytes):
            line = line.decode()
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                seqs[name] = seq.upper()
            name = l[1:]
            seq = ''
        else:
            seq += l
    if name is not None and len(seq) > 0:
        seqs[name] = seq.upper()
    return seqs


def parse_header(header):
    parts = header.split('|')
    return {
        'name': parts[0] if len(parts) > 0 else '',
        'accession': parts[1] if len(parts) > 1 else '',
        'country': parts[2] if len(parts) > 2 else '',
        'admin1': parts[3] if len(parts) > 3 else '',
        'admin2': parts[4] if len(parts) > 4 else '',
        'date': parts[5] if len(parts) > 5 else '',
    }


def compute_metrics(seq):
    length = len(seq)
    n_count = seq.count('N')
    gap_count = seq.count('-')
    acgt = sum(1 for c in seq if c in 'ACGT')
    called = acgt
    # Effective length: first called base to last called base
    first_called = None
    last_called = None
    for i, c in enumerate(seq):
        if c in 'ACGT':
            if first_called is None:
                first_called = i
            last_called = i
    if first_called is not None:
        effective_length = last_called - first_called + 1
    else:
        effective_length = 0
    return {
        'length': length,
        'effective_length': effective_length,
        'n_count': n_count,
        'n_fraction': n_count / length if length > 0 else 0,
        'gap_count': gap_count,
        'acgt_count': acgt,
        'called_fraction': called / length if length > 0 else 0,
    }


def check_date(date_str):
    if re.match(r'^\d{4}-\d{2}-\d{2}$', date_str):
        return 'complete'
    if re.match(r'^\d{4}-\d{2}$', date_str):
        return 'year-month'
    if re.match(r'^\d{4}$', date_str):
        return 'year-only'
    return 'missing'


def load_original_dates(metadata_path):
    """Load original sampleCollectionDate from the Pathoplexus metadata TSV."""
    dates = {}
    with open(metadata_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            acc = row.get('accessionVersion', '')
            date = row.get('sampleCollectionDate', '')
            dates[acc] = date
    return dates


def main():
    parser = argparse.ArgumentParser(description="Basic genome QC report")
    parser.add_argument("-a", "--alignment", required=True, help="Cleaned alignment FASTA")
    parser.add_argument("-m", "--metadata", required=True, help="Pathoplexus metadata TSV")
    parser.add_argument("-o", "--output", required=True, help="Output report file path")
    parser.add_argument("--max-n-fraction", type=float, default=0.10,
                        help="Flag genomes with N fraction above this (default: 0.10)")
    parser.add_argument("--min-completeness", type=float, default=0.80,
                        help="Flag genomes with effective length < this fraction of reference (default: 0.80)")
    parser.add_argument("--min-called-fraction", type=float, default=0.50,
                        help="Flag genomes with called bases < this fraction (default: 0.50)")
    args = parser.parse_args()

    seqs = readFASTA(args.alignment)
    original_dates = load_original_dates(args.metadata)

    records = []
    for header, seq in seqs.items():
        info = parse_header(header)
        metrics = compute_metrics(seq)
        original_date = original_dates.get(info['accession'], info['date'])
        date_status = check_date(original_date)
        records.append({**info, **metrics, 'date_status': date_status, 'original_date': original_date})

    eff_lengths = [r['effective_length'] for r in records if r['effective_length'] > 0]
    median_eff_len = statistics.median(eff_lengths) if eff_lengths else 0
    min_eff_len = args.min_completeness * median_eff_len

    for r in records:
        flags = []
        if r['n_fraction'] > args.max_n_fraction:
            flags.append('high_Ns')
        if r['effective_length'] < min_eff_len:
            flags.append('short')
        if r['called_fraction'] < args.min_called_fraction:
            flags.append('low_called')
        if r['date_status'] != 'complete':
            flags.append('incomplete_date')
        r['flags'] = flags

    flagged = [r for r in records if r['flags']]
    flagged.sort(key=lambda r: (-r['n_fraction'], r['name']))

    counts = {}
    for r in flagged:
        for f in r['flags']:
            counts[f] = counts.get(f, 0) + 1

    hdr = f"{'name':<16} {'accession':<16} {'location':<22} {'date':<12} {'eff_len':>7} {'called%':>7} {'N%':>6}  flags"
    sep = '-' * len(hdr)

    lines = []
    lines.append("BASIC GENOME QC REPORT")
    lines.append("=" * 50)
    lines.append(f"Alignment: {len(records)} sequences, {median_eff_len:.0f} bp median effective length")
    lines.append(f"Thresholds: N fraction > {args.max_n_fraction:.0%}, effective length < {args.min_completeness:.0%} of median ({min_eff_len:.0f} bp), called < {args.min_called_fraction:.0%}")
    lines.append("")

    if flagged:
        lines.append(f"FLAGGED GENOMES ({len(flagged)} of {len(records)}):")
        lines.append(hdr)
        lines.append(sep)
        for r in flagged:
            loc = f"{r['admin1']}/{r['admin2']}" if r['admin2'] else r['admin1']
            lines.append(
                f"{r['name']:<16} {r['accession']:<16} {loc:<22} {r['original_date']:<12} "
                f"{r['effective_length']:>7} {r['called_fraction']:>6.1%} {r['n_fraction']:>5.1%}  {', '.join(r['flags'])}"
            )
        lines.append("")
    else:
        lines.append("No genomes flagged.")
        lines.append("")

    lines.append("SUMMARY:")
    for flag_name in ['high_Ns', 'short', 'low_called', 'incomplete_date']:
        c = counts.get(flag_name, 0)
        lines.append(f"  {flag_name:<20} {c} genomes")
    lines.append(f"  {'total flagged':<20} {len(flagged)} of {len(records)}")

    report = '\n'.join(lines)

    with open(args.output, 'w') as f:
        f.write(report + '\n')

    print(report)


if __name__ == '__main__':
    main()
