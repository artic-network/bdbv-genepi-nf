#!/usr/bin/env python3
"""List new duplicate_sample flags from health checks that are not yet in the pending exclusion list."""
import argparse
from pathlib import Path

from health_checks import run_all_checks


def read_exclude_entries(path):
    entries = set()
    if not path or not Path(path).exists():
        return entries
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                entries.add(line.split()[0])
    return entries


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--current', required=True, help='Current pipeline output directory')
parser.add_argument('--previous', default=None, help='Previous pipeline output directory')
parser.add_argument('--prefix', required=True, help='File prefix pattern')
parser.add_argument('--full', action='store_true', help='Check all genomes, not just new ones')
parser.add_argument('--exclude-pending', default=None, help='Existing pending exclusion list')
parser.add_argument('-o', '--output', required=True, help='Output file of sample IDs to add')
args = parser.parse_args()

existing = read_exclude_entries(args.exclude_pending)
results = run_all_checks(args.current, args.previous, args.prefix, full=args.full)
new_pending = []
for entry in results.get('flagged_sequences', []):
    flag_names = {f['flag'] for f in entry['flags']}
    if 'duplicate_sample' in flag_names and entry['sample_id'] not in existing and entry['sample_id'] not in new_pending:
        new_pending.append(entry['sample_id'])

with open(args.output, 'w') as f:
    f.write('# duplicate_sample flags to add to data/exclude_pending.txt\n')
    for sid in new_pending:
        f.write(f'{sid}\n')
print(f"{len(new_pending)} new pending exclusions written to {args.output}")
