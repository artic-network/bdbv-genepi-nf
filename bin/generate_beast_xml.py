#!/usr/bin/env python3
"""Generate BEAST XML(s) from a template by swapping in taxa, alignment, and SkyGrid params."""

import argparse
import math
import re
import sys
from datetime import datetime
from pathlib import Path


def parse_date(s):
    for fmt in ('%Y-%m-%d', '%Y-%m', '%Y'):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    return None


def date_to_fracyear(d):
    year = d.year
    jan1 = datetime(year, 1, 1).date()
    days_in_year = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    return year + (d - jan1).days / days_in_year


def header_to_date_str(header):
    parts = header.split('|')
    return parts[5] if len(parts) >= 6 else parts[-1]


def header_to_fracyear(header):
    d = parse_date(header_to_date_str(header))
    return date_to_fracyear(d) if d else None


def parse_alignment(fasta_path):
    seqs = {}
    order = []
    current = None
    chunks = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current is not None:
                    seqs[current] = ''.join(chunks)
                current = line[1:]
                order.append(current)
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        seqs[current] = ''.join(chunks)
    return seqs, order


def compute_skygrid_params(headers, anchor_date):
    """Compute SkyGrid cutOff and numGridPoints for weekly intervals."""
    anchor = parse_date(anchor_date)
    if anchor is None:
        print(f"ERROR: could not parse anchor date '{anchor_date}'", file=sys.stderr)
        sys.exit(1)

    latest = None
    for h in headers:
        d = parse_date(header_to_date_str(h))
        if d and (latest is None or d > latest):
            latest = d

    week_years = 7.0 / 365.25
    weeks = round((latest - anchor).days / 7)
    if weeks < 1:
        print(f"ERROR: anchor date {anchor_date} is after the latest sample {latest}",
              file=sys.stderr)
        sys.exit(1)
    cutoff = weeks * week_years
    num_grid_points = weeks - 1

    print(f"  SkyGrid: anchor={anchor_date}, latest={latest}, "
          f"weeks={weeks}, numGridPoints={num_grid_points}, cutOff={cutoff:.10f}")

    return {'numGridPoints': num_grid_points, 'cutOff': cutoff}


def generate_xml(template_path, headers, seqs, skygrid_params,
                 new_prefix, chain_length, log_every, output_path):
    """Generate a BEAST XML from a template by swapping taxa, alignment, and params."""
    text = open(template_path).read()
    n = len(headers)

    m = re.search(r'<log\s+id="fileLog"[^>]*fileName="([^"]+)\.log"', text)
    if not m:
        print(f'  WARNING: could not find fileLog in {template_path}, skipping')
        return
    old_prefix = m.group(1)
    text = text.replace(old_prefix, new_prefix)

    taxa_lines = ['\t<taxa id="taxa">']
    for h in headers:
        fy = header_to_fracyear(h)
        if fy is None:
            continue
        taxa_lines.append(f'\t\t<taxon id="{h}">')
        taxa_lines.append(f'\t\t\t<date value="{fy}" direction="forwards" units="years"/>')
        taxa_lines.append(f'\t\t</taxon>')
    taxa_lines.append('\t</taxa>')
    taxa_block = '\n'.join(taxa_lines)
    text = re.sub(r'\t<taxa id="taxa">.*?</taxa>', taxa_block, text, flags=re.DOTALL)

    aln_lines = ['\t<alignment id="alignment" dataType="nucleotide">']
    for h in headers:
        aln_lines.append(f'\t\t<sequence>')
        aln_lines.append(f'\t\t\t<taxon idref="{h}"/>')
        aln_lines.append(f'\t\t\t{seqs[h]}')
        aln_lines.append(f'\t\t</sequence>')
    aln_lines.append('\t</alignment>')
    aln_block = '\n'.join(aln_lines)
    text = re.sub(r'\t<alignment id="alignment"[^>]*>.*?</alignment>', aln_block, text, flags=re.DOTALL)

    text = re.sub(r'<!-- ntax=\d+', f'<!-- ntax={n}', text)

    if skygrid_params and 'skygrid.numGridPoints' in text:
        ngp = skygrid_params['numGridPoints']
        cutoff = skygrid_params['cutOff']
        text = re.sub(
            r'(<parameter\s+id="skygrid\.logPopSize"\s+dimension=")\d+(")',
            rf'\g<1>{ngp + 1}\2', text)
        text = re.sub(
            r'(<parameter\s+id="skygrid\.numGridPoints"\s+value=")[^"]+(")',
            rf'\g<1>{ngp}.0\2', text)
        text = re.sub(
            r'(<parameter\s+id="skygrid\.cutOff"\s+value=")[^"]+(")',
            rf'\g<1>{cutoff:.10f}\2', text)

    text = re.sub(
        r'(<subtreeLeap\s+size="[^"]*"\s+weight=")\S+(")',
        rf'\g<1>{n}\2', text)
    text = re.sub(
        r'(<fixedHeightSubtreePruneRegraft\s+weight=")\S+(")',
        rf'\g<1>{round(n / 10, 1)}\2', text)
    text = re.sub(
        r'(<mutationHistoryGibbsOperator\s+weight=")\S+(")',
        rf'\g<1>{round(n / 10, 1)}\2', text)

    text = re.sub(
        r'(<mcmc\s+id="mcmc"\s+chainLength=")\d+(")',
        rf'\g<1>{chain_length}\2', text)

    text = re.sub(
        r'(<log\s+id="fileLog"\s+logEvery=")\d+(")',
        rf'\g<1>{log_every}\2', text)
    text = re.sub(
        r'(<logTree\s+id="treeFileLog"\s+logEvery=")\d+(")',
        rf'\g<1>{log_every}\2', text)

    burnin = chain_length // 10
    text = re.sub(r'(burnin=")\d+(")', rf'\g<1>{burnin}\2', text)

    text = re.sub(
        r'(checkpointEvery=")\d+(")',
        rf'\g<1>{chain_length // 10}\2', text)
    text = re.sub(
        r'(checkpointFinal=")\d+(")',
        rf'\g<1>{chain_length}\2', text)

    with open(output_path, 'w') as f:
        f.write(text)
    print(f'  Wrote {output_path} ({n} taxa)')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-a', '--alignment', required=True,
                        help='Input FASTA alignment (pipeline-formatted headers)')
    parser.add_argument('-o', '--outdir', required=True, type=Path,
                        help='Output directory for generated XMLs')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Prefix for output file names and BEAST log files')
    parser.add_argument('--template-egc', default=None,
                        help='EGC (exponential growth coalescent) template XML')
    parser.add_argument('--template-sg', default=None,
                        help='SkyGrid template XML')
    parser.add_argument('--template-lgc', default=None,
                        help='LGC (logistic growth coalescent) template XML')
    parser.add_argument('--template-mat', default=None,
                        help='MAT (mutation-annotated tree) SkyGrid template XML')
    parser.add_argument('--anchor-date', default=None,
                        help='SkyGrid anchor date (YYYY-MM-DD). Required for SkyGrid templates.')
    parser.add_argument('--chain-length', type=int, default=500_000_000,
                        help='MCMC chain length (default: 500M)')
    parser.add_argument('--log-every', type=int, default=50_000,
                        help='Log sampling interval (default: 50000)')
    args = parser.parse_args()

    if (args.template_sg or args.template_mat) and not args.anchor_date:
        print("ERROR: --anchor-date is required when using --template-sg or --template-mat", file=sys.stderr)
        sys.exit(1)

    if not args.template_egc and not args.template_sg and not args.template_lgc and not args.template_mat:
        print("ERROR: provide at least one of --template-egc, --template-sg, or --template-lgc", file=sys.stderr)
        sys.exit(1)

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading alignment: {args.alignment}")
    seqs, headers = parse_alignment(args.alignment)
    print(f"  {len(headers)} sequences, {len(seqs[headers[0]])} sites")

    if args.template_egc:
        print(f"\nGenerating EGC XML from {args.template_egc}")
        out_path = args.outdir / f"{args.prefix}.GTR_EGC.xml"
        generate_xml(Path(args.template_egc), headers, seqs,
                     None, args.prefix, args.chain_length, args.log_every, out_path)

    if args.template_sg:
        print(f"\nGenerating SkyGrid XML from {args.template_sg}")
        sg_params = compute_skygrid_params(headers, args.anchor_date)
        out_path = args.outdir / f"{args.prefix}.GTR_SG.xml"
        generate_xml(Path(args.template_sg), headers, seqs,
                     sg_params, args.prefix, args.chain_length, args.log_every, out_path)

    if args.template_lgc:
        print(f"\nGenerating LGC XML from {args.template_lgc}")
        out_path = args.outdir / f"{args.prefix}.GTR_LGC.xml"
        generate_xml(Path(args.template_lgc), headers, seqs,
                     None, args.prefix, args.chain_length, args.log_every, out_path)

    if args.template_mat:
        print(f"\nGenerating MAT SkyGrid XML from {args.template_mat}")
        sg_params = compute_skygrid_params(headers, args.anchor_date)
        out_path = args.outdir / f"{args.prefix}.GTR_SG.mat.xml"
        generate_xml(Path(args.template_mat), headers, seqs,
                     sg_params, args.prefix, args.chain_length, args.log_every, out_path)

    print("\nDone.")


if __name__ == '__main__':
    main()
