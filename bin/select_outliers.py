#!/usr/bin/env python3
"""Keep only temporal outliers with |z| >= threshold from a temporal_root.py outliers file."""
import argparse
import re

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--input", required=True, help="Detected outliers file (<prefix>.outliers.txt)")
parser.add_argument("-t", "--threshold", type=float, required=True, help="Minimum |z| to remove")
parser.add_argument("-o", "--output", required=True, help="Output outlier exclusion file")
args = parser.parse_args()

kept = []
total = 0
with open(args.input) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        total += 1
        m = re.search(r'z=[+-]?([\d.]+)', line)
        if m and float(m.group(1)) >= args.threshold:
            kept.append(line)

with open(args.output, 'w') as f:
    f.write(f"# Outliers auto-removed at |z| >= {args.threshold}\n")
    f.writelines(kept)
print(f"Auto-removal: {len(kept)} outliers with |z| >= {args.threshold} (from {total} detected)")
