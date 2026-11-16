#!/usr/bin/env python3
"""Check that a FASTA file is an alignment (all sequences the same length).

Exits with status 1 and a clear message if not, so a bad input fails at the start of the
pipeline instead of when IQ-TREE reads it.
"""
import argparse
import sys
from collections import Counter

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("fasta", help="FASTA file to check")
parser.add_argument("--mafft-hint", action="store_true",
                    help="Suggest running with --mafft if the sequences are not aligned")
args = parser.parse_args()

lengths = {}
name = None
with open(args.fasta) as f:
    for line in f:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            name = line[1:].strip()
            if name in lengths:
                print(f"ERROR: duplicate sequence name in {args.fasta}: {name}", file=sys.stderr)
                sys.exit(1)
            lengths[name] = 0
        elif name is not None:
            lengths[name] += len(line.strip())

if not lengths:
    print(f"ERROR: no sequences found in {args.fasta}", file=sys.stderr)
    sys.exit(1)

counts = Counter(lengths.values())
if len(counts) > 1:
    common_len, common_n = counts.most_common(1)[0]
    odd = [n for n, l in lengths.items() if l != common_len]
    print(f"ERROR: {args.fasta} is not an alignment: sequences have {len(counts)} different lengths "
          f"(shortest {min(counts)}, longest {max(counts)}).", file=sys.stderr)
    print(f"  {common_n} of {len(lengths)} sequences are {common_len} long; "
          f"{len(odd)} differ, e.g. {', '.join(odd[:3])}", file=sys.stderr)
    if args.mafft_hint:
        print("  The input looks unaligned. Rerun with --mafft to align it with MAFFT first.", file=sys.stderr)
    sys.exit(1)

print(f"Alignment check passed: {len(lengths)} sequences, all {next(iter(counts))} long")
