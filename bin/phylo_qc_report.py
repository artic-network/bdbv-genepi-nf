#!/usr/bin/env python3
import argparse
import csv
import re


GP_START = 5875
GP_END = 8279


def parse_header_name(child_field):
    """Extract the short name (e.g. BIA-1299) from a branch_snps child field like
    'BIA-1299|PP_00764QW.1|DRC|Ituri|Rwampara|2026-06-04'."""
    return child_field.split('|')[0]


def parse_header_full(header):
    parts = header.split('|')
    return {
        'name': parts[0] if len(parts) > 0 else '',
        'accession': parts[1] if len(parts) > 1 else '',
        'admin1': parts[3] if len(parts) > 3 else '',
        'admin2': parts[4] if len(parts) > 4 else '',
        'date': parts[5] if len(parts) > 5 else '',
    }


def is_tip(child_field):
    return '|' in child_field


def parse_branch_snps(path):
    """Returns {tip_label: [(site, snp_type), ...]} for terminal branches only."""
    tips = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            child = row['child']
            if not is_tip(child):
                continue
            name = parse_header_name(child)
            site = int(row['site'])
            snp = row['snp']
            tips.setdefault(name, {'label': child, 'snps': []})
            tips[name]['snps'].append((site, snp))
    return tips


def parse_phylo_flags(path):
    """Returns set of tip names flagged as long_branch."""
    long_branch = set()
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['mutation_type'] == 'long_branch':
                for tip in row['present_in'].split(';'):
                    tip = tip.strip()
                    if tip:
                        name = tip.split('|')[0]
                        long_branch.add(name)
    return long_branch


def count_snp_types(snps):
    tc = sum(1 for _, s in snps if s == 'T->C')
    ag = sum(1 for _, s in snps if s == 'A->G')
    other = len(snps) - tc - ag
    gp_tc = sum(1 for site, s in snps if s == 'T->C' and GP_START <= site <= GP_END)
    nongp_tc = tc - gp_tc
    return {
        'total': len(snps),
        'tc': tc,
        'ag': ag,
        'other': other,
        'gp_tc': gp_tc,
        'nongp_tc': nongp_tc,
    }


def parse_basic_qc(path):
    """Parse basic QC report to get per-genome N fraction. Returns {name: n_fraction}."""
    if path is None:
        return {}
    n_fracs = {}
    try:
        with open(path) as f:
            in_table = False
            for line in f:
                line = line.rstrip()
                if line.startswith('name') and 'N%' in line:
                    in_table = True
                    continue
                if in_table and line.startswith('-'):
                    continue
                if in_table and line.strip() == '':
                    in_table = False
                    continue
                if in_table:
                    parts = line.split()
                    if len(parts) >= 7:
                        name = parts[0]
                        n_pct_str = parts[6]
                        try:
                            n_fracs[name] = float(n_pct_str.rstrip('%')) / 100
                        except ValueError:
                            pass
    except FileNotFoundError:
        pass
    return n_fracs


def load_metadata_ns(path):
    """Load totalUnknownNucs from metadata TSV. Returns {submissionId: n_count}."""
    ns = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sid = row.get('submissionId', '')
            try:
                n_count = int(row.get('totalUnknownNucs', 0))
            except (ValueError, TypeError):
                n_count = 0
            if sid:
                ns[sid] = n_count
    return ns


def assign_concern(snp_counts, is_long_branch, n_fraction, tc_threshold, snp_threshold):
    tc = snp_counts['tc']
    ag = snp_counts['ag']
    other = snp_counts['other']
    total = snp_counts['total']
    high_ns = n_fraction is not None and n_fraction > 0.10

    reasons = []

    if is_long_branch:
        reasons.append('long_branch')
    if tc > 0:
        region = 'GP-only' if snp_counts['nongp_tc'] == 0 else 'genome-wide'
        reasons.append(f'{tc} T>C ({region})')
    if ag > 0:
        reasons.append(f'{ag} A>G')
    if other > 0:
        reasons.append(f'{other} other')
    if high_ns:
        reasons.append('high_Ns')

    # examine: most suspicious
    if is_long_branch and tc > tc_threshold:
        return 'examine', ' + '.join(reasons)
    if tc > tc_threshold and other > 0:
        return 'examine', ' + '.join(reasons)
    if high_ns and tc > tc_threshold:
        return 'examine', ' + '.join(reasons)
    if total > snp_threshold and tc > tc_threshold:
        return 'examine', ' + '.join(reasons)
    if tc > tc_threshold:
        return 'examine', ' + '.join(reasons)

    # watch: worth noting
    if is_long_branch:
        return 'watch', ' + '.join(reasons)
    if 0 < tc <= tc_threshold:
        return 'watch', ' + '.join(reasons)
    if high_ns:
        return 'watch', ' + '.join(reasons)
    if ag > 2:
        return 'watch', ' + '.join(reasons)
    if total > snp_threshold:
        return 'watch', ' + '.join(reasons)

    return 'ok', ''


def main():
    parser = argparse.ArgumentParser(description="Phylogenetic QC report")
    parser.add_argument("-f", "--phylo-flags", required=True, help="Raccoon phylo_flags.csv")
    parser.add_argument("-b", "--branch-snps", required=True, help="Raccoon branch_snps.reconstruction.csv")
    parser.add_argument("-m", "--metadata", required=True, help="Pathoplexus metadata TSV")
    parser.add_argument("--basic-qc", default=None, help="Basic QC report for N content cross-reference")
    parser.add_argument("-o", "--output", required=True, help="Output report file path")
    parser.add_argument("--tc-threshold", type=int, default=2,
                        help="Flag genomes with T>C count above this (default: 2)")
    parser.add_argument("--snp-threshold", type=int, default=5,
                        help="Flag terminal branches with total SNPs above this (default: 5)")
    args = parser.parse_args()

    long_branch_tips = parse_phylo_flags(args.phylo_flags)
    tip_snps = parse_branch_snps(args.branch_snps)
    basic_ns = parse_basic_qc(args.basic_qc)
    meta_ns = load_metadata_ns(args.metadata)

    all_tips = set()
    tip_info = {}

    for name, data in tip_snps.items():
        all_tips.add(name)
        info = parse_header_full(data['label'])
        tip_info[name] = info

    records = []
    for name in sorted(all_tips):
        data = tip_snps.get(name, {'snps': []})
        snps = data.get('snps', [])
        snp_counts = count_snp_types(snps)
        is_lb = name in long_branch_tips
        n_frac = basic_ns.get(name)
        meta_n = meta_ns.get(name, 0)
        info = tip_info.get(name, {})

        concern, reason = assign_concern(snp_counts, is_lb, n_frac, args.tc_threshold, args.snp_threshold)

        records.append({
            'name': name,
            'accession': info.get('accession', ''),
            'admin1': info.get('admin1', ''),
            'admin2': info.get('admin2', ''),
            'date': info.get('date', ''),
            'meta_ns': meta_n,
            'n_frac': n_frac,
            'total': snp_counts['total'],
            'tc': snp_counts['tc'],
            'ag': snp_counts['ag'],
            'other': snp_counts['other'],
            'gp_tc': snp_counts['gp_tc'],
            'nongp_tc': snp_counts['nongp_tc'],
            'long_branch': is_lb,
            'concern': concern,
            'reason': reason,
        })

    concern_order = {'examine': 0, 'watch': 1, 'ok': 2}
    records.sort(key=lambda r: (concern_order[r['concern']], -r['total'], r['name']))

    flagged = [r for r in records if r['concern'] != 'ok']

    n_examine = sum(1 for r in records if r['concern'] == 'examine')
    n_watch = sum(1 for r in records if r['concern'] == 'watch')
    n_ok = sum(1 for r in records if r['concern'] == 'ok')

    hdr = (f"{'name':<16} {'location':<22} {'date':<12} {'Ns':>6} "
           f"{'SNPs':>5} {'T>C':>4} {'A>G':>4} {'othr':>4} {'LB':>3}  {'concern':<8} reason")
    sep = '-' * len(hdr)

    lines = []
    lines.append("PHYLOGENETIC QC REPORT")
    lines.append("=" * 50)
    lines.append(f"Genomes with terminal branch SNPs: {len(all_tips)}")
    lines.append(f"Long branch flags: {len(long_branch_tips)}")
    lines.append(f"Thresholds: T>C > {args.tc_threshold}, total SNPs > {args.snp_threshold}")
    lines.append("")

    if flagged:
        lines.append(f"FLAGGED GENOMES ({len(flagged)}):")
        lines.append(hdr)
        lines.append(sep)
        for r in flagged:
            loc = f"{r['admin1']}/{r['admin2']}" if r['admin2'] else r['admin1']
            lb = 'yes' if r['long_branch'] else ''
            ns_str = f"{r['meta_ns']:>6}" if r['meta_ns'] else '     .'
            lines.append(
                f"{r['name']:<16} {loc:<22} {r['date']:<12} {ns_str} "
                f"{r['total']:>5} {r['tc']:>4} {r['ag']:>4} {r['other']:>4} {lb:>3}  "
                f"{r['concern']:<8} {r['reason']}"
            )
        lines.append("")
    else:
        lines.append("No genomes flagged.")
        lines.append("")

    lines.append("SUMMARY:")
    lines.append(f"  examine:  {n_examine} genomes")
    lines.append(f"  watch:    {n_watch} genomes")
    lines.append(f"  ok:       {n_ok} genomes")

    report = '\n'.join(lines)

    with open(args.output, 'w') as f:
        f.write(report + '\n')

    print(report)


if __name__ == '__main__':
    main()
