#!/usr/bin/env python3
"""Check that a FASTA file is an alignment (all sequences the same length).

Exits with status 1 and a clear message if not, so a bad input fails at the start of the
pipeline instead of when IQ-TREE reads it. With --quiet nothing is printed and only the exit
status is used (0 = aligned, 1 = not), so a script can ask "is this already aligned?".
"""
import argparse
import sys
from collections import Counter

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("fasta", help="FASTA file to check")
parser.add_argument("--mafft-hint", action="store_true",
                    help="Suggest running with --mafft if the sequences are not aligned")
parser.add_argument("--quiet", action="store_true",
                    help="Print nothing; exit 0 if the file is an alignment and 1 if not")
args = parser.parse_args()


def say(message, stream=sys.stderr):
    if not args.quiet:
        print(message, file=stream)

lengths = {}
name = None
with open(args.fasta) as f:
    for line in f:
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            name = line[1:].strip()
            if name in lengths:
                say(f"ERROR: duplicate sequence name in {args.fasta}: {name}")
                sys.exit(1)
            lengths[name] = 0
        elif name is not None:
            lengths[name] += len(line.strip())

if not lengths:
    say(f"ERROR: no sequences found in {args.fasta}")
    sys.exit(1)

counts = Counter(lengths.values())
if len(counts) > 1:
    common_len, common_n = counts.most_common(1)[0]
    odd = [n for n, l in lengths.items() if l != common_len]
    say(f"ERROR: {args.fasta} is not an alignment: sequences have {len(counts)} different lengths "
        f"(shortest {min(counts)}, longest {max(counts)}).")
    say(f"  {common_n} of {len(lengths)} sequences are {common_len} long; "
        f"{len(odd)} differ, e.g. {', '.join(odd[:3])}")
    if args.mafft_hint:
        say("  The input looks unaligned. Rerun with --mafft to align it with MAFFT first.")
    sys.exit(1)

say(f"Alignment check passed: {len(lengths)} sequences, all {next(iter(counts))} long", sys.stdout)
