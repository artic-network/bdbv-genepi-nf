import pandas as pd
from datetime import datetime
from pathlib import Path
import argparse
from gzip import open as gopen

# read FASTA stream and return (ID,seq) dictionary
def readFASTA(path):
    if path.endswith('.gz'):
        stream = gopen(path)
    else:
        stream = open(path)
    seqs = {}
    name = None
    seq = ''
    for line in stream:
        if isinstance(line,bytes):
            line = line.decode()
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                assert len(seq) != 0, "Malformed FASTA"
                seqs[name] = seq.upper()
            name = l[1:]
            assert name not in seqs, "Duplicate sequence ID: %s" % name
            seq = ''
        else:
            seq += l
    assert name is not None and len(seq) != 0, "Malformed FASTA"
    seqs[name] = seq.upper()
    return seqs


def writeFASTA(seqs, path):
    with open(path, 'w') as f:
        for k in seqs:
            f.write('>%s\n' % k)
            f.write('%s\n' % seqs[k])
            

country_aliases = {
    'Democratic Republic of the Congo': 'DRC',
    'Democratic Republic of Congo': 'DRC',
    'Democratic Republic of the Congo (DRC)': 'DRC',
    'democratic republic of the congo': 'DRC',
}

# argparse
parser = argparse.ArgumentParser(description="Clean up Pathoplexus data")
parser.add_argument("-f", "--fasta", type=str, required=True, help="Path to the Pathoplexus fasta file.")
parser.add_argument("-m", "--metadata", type=str, required=True, help="Path to the Pathoplexus metadata file.")
parser.add_argument("-a", "--aliases_file", type=str, required=False, default="data/aliases.csv", help="Path to the aliases file containing mappings for DRC provinces and health zones.")
parser.add_argument("-ts", "--trim_start", type=int, required=False, help="First position to keep (1-based, inclusive). Positions before this are masked to N.")
parser.add_argument("-te", "--trim_end", type=int, required=False, help="Last position to keep (1-based, inclusive). Positions after this are masked to N.")
parser.add_argument("--mask_sites", type=str, required=False, help="Comma-separated 1-based positions to mask to N (e.g. 5246,5249,5251).")
parser.add_argument("--mask_sites_file", type=str, required=False, help="File with 1-based positions to mask, one per line (# comments allowed).")
parser.add_argument("-e", "--exclude-files", type=str, nargs='*', default=[], help="One or more text files with accession IDs to exclude (one per line, # comments allowed).")
parser.add_argument("-p", "--prefix", type=str, required=True, help="Prefix for the output file.")
parser.add_argument("--min-date", type=str, default=None, help="Remove sequences with date before this (YYYY-MM-DD). Sequences in the same year are flagged with a warning.")

args = parser.parse_args()

# Parse mask sites (1-based → 0-based indices)
mask_positions = set()
if args.mask_sites:
    for s in args.mask_sites.split(','):
        mask_positions.add(int(s.strip()) - 1)
if args.mask_sites_file:
    with open(args.mask_sites_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            mask_positions.add(int(line.split()[0]) - 1)
if mask_positions:
    print(f"Masking {len(mask_positions)} sites to N: {','.join(str(p+1) for p in sorted(mask_positions))}")

# Read inputs
fasta = readFASTA(args.fasta)
meta = pd.read_csv(args.metadata, sep="\t")
aliases = pd.read_csv(args.aliases_file)

# Build admin2 alias lookup from aliases file (case-insensitive)
admin2_aliases = dict(zip(aliases['observed_name'].str.lower(), aliases['canonical_nom']))

# Standardize metadata fields
meta['geoLocCountry_clean'] = meta['geoLocCountry'].map(country_aliases).fillna(meta['geoLocCountry'])
meta['geoLocAdmin1_clean'] = meta['geoLocAdmin1'].fillna('')
meta['geoLocAdmin2_clean'] = meta['geoLocAdmin2'].fillna('').apply(
    lambda x: admin2_aliases.get(x.lower(), x) if x else '')

# Deduplicate: keep only the highest version per accession
meta['_acc_base'] = meta['accessionVersion'].str.rsplit('.', n=1).str[0]
meta['_acc_ver'] = meta['accessionVersion'].str.rsplit('.', n=1).str[1].astype(int)
meta = meta.sort_values('_acc_ver').drop_duplicates(subset='_acc_base', keep='last')
meta = meta.drop(columns=['_acc_base', '_acc_ver'])

# Deduplicate FASTA: keep only the highest version per accession
fasta_by_base = {}
for acc_id in fasta:
    base = acc_id.rsplit('.', 1)[0]
    ver = int(acc_id.rsplit('.', 1)[1])
    if base not in fasta_by_base or ver > fasta_by_base[base][0]:
        fasta_by_base[base] = (ver, acc_id)
kept_accs = {v[1] for v in fasta_by_base.values()}
n_deduped = len(fasta) - len(kept_accs)
if n_deduped > 0:
    print(f"Version dedup: {len(fasta)} → {len(kept_accs)} sequences ({n_deduped} older versions removed)")
fasta = {k: v for k, v in fasta.items() if k in kept_accs}

# Index metadata by accessionVersion for FASTA lookup
meta_by_acc = meta.set_index('accessionVersion')

# Build new FASTA with standardized headers
# Format: submissionId|accessionVersion|geoLocCountry|geoLocAdmin1|geoLocAdmin2|collectionDate
out_seqs = {}
missing_meta = []
unknown_admin2 = set()

# PP_006XXY5: withdrawn sequence on Pathoplexus (replaced by corrected submission)
to_remove = set(['PP_006XXY5'])

for ef in args.exclude_files:
    with open(ef) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            token = line.split()[0]
            if '|' in token:
                acc = token.split('|')[1].split('.')[0]
            else:
                acc = token.split('.')[0]
            to_remove.add(acc)

print(f"Excluding {len(to_remove)} accessions ({len(to_remove) - 1} from exclude files + 1 default)")

for acc_id, seq in fasta.items():
    if acc_id.split('.')[0] in to_remove:
        continue
    
    if acc_id not in meta_by_acc.index:
        missing_meta.append(acc_id)
        continue

    row = meta_by_acc.loc[acc_id]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]

    submission_id = str(row.get('submissionId', ''))
    country = str(row['geoLocCountry_clean'])
    admin1 = str(row['geoLocAdmin1_clean'])
    admin2 = str(row['geoLocAdmin2_clean'])
    date = str(row.get('sampleCollectionDate', ''))
    
    if len(date) < 10:
        print(f"WARNING: Incomplete date for {acc_id}: {date}")
        print(f"  Using sampleReceivedDate instead: {row.get('sampleReceivedDate', '')}")
        date = str(row.get('sampleReceivedDate', ''))

    if admin2 and admin2.lower() not in admin2_aliases and admin2 not in admin2_aliases.values() and admin2 not in ('', 'nan'):
        unknown_admin2.add(admin2)

    seq = seq.upper()
    start = 1
    end = len(seq)
    if args.trim_start is not None:
        start = args.trim_start
    if args.trim_end is not None:
        end = args.trim_end
    if start > 1 or end < len(seq):
        seq = 'N' * (start - 1) + seq[start - 1:end] + 'N' * (len(seq) - end)

    for idx in mask_positions:
        if idx < len(seq):
            seq = seq[:idx] + 'N' + seq[idx+1:]

    # replace gaps with N's
    seq = seq.replace('-', 'N')
    
    # replace spaces in admin1/admin2 with hyphens
    admin1 = admin1.replace(' ', '-') if admin1 else ''
    admin2 = admin2.replace(' ', '-') if admin2 else ''
        
    new_header = f"{submission_id}|{acc_id}|{country}|{admin1}|{admin2}|{date}"
    out_seqs[new_header] = seq

# Min-date filtering
if args.min_date:
    min_dt = datetime.strptime(args.min_date, '%Y-%m-%d')
    min_year = min_dt.year
    flagged = []
    old_removed = []
    bad_date_removed = []
    to_keep = {}
    for header, seq in out_seqs.items():
        date_str = header.split('|')[-1]
        if date_str in ('', 'nan') or len(date_str) < 10:
            bad_date_removed.append(header)
            continue
        try:
            seq_dt = datetime.strptime(date_str[:10], '%Y-%m-%d')
        except ValueError:
            bad_date_removed.append(header)
            continue
        if seq_dt < min_dt:
            if seq_dt.year == min_year:
                flagged.append((header, date_str))
            else:
                old_removed.append(header)
            continue
        to_keep[header] = seq

    if old_removed:
        print(f"\nDate filter: removed {len(old_removed)} sequences from before {min_year}")
    if bad_date_removed:
        print(f"Date filter: removed {len(bad_date_removed)} sequences with unparseable/missing dates")
    if flagged:
        print(f"\nWARNING: {len(flagged)} sequences dated {min_year} but before {args.min_date} (removed):")
        for header, date_str in flagged:
            fields = header.split('|')
            print(f"  {fields[0]} ({fields[1]}) {fields[2]} {date_str}")
        flagged_path = f"{args.prefix}.flagged_early_dates.txt"
        with open(flagged_path, 'w') as f:
            f.write(f"# Sequences dated {min_year} but before {args.min_date}\n")
            for header, date_str in flagged:
                f.write(f"{header}\n")
        print(f"  Written to {flagged_path}")

    n_before = len(out_seqs)
    out_seqs = to_keep
    n_removed = n_before - len(out_seqs)
    if n_removed > 0:
        print(f"Date filter: {n_before} → {len(out_seqs)} sequences ({n_removed} removed)")

# Report
print(f"Input: {len(fasta)} sequences, {len(meta)} metadata rows")
print(f"Output: {len(out_seqs)} sequences with standardized headers")

if missing_meta:
    print(f"\nWARNING: {len(missing_meta)} sequences not found in metadata:")
    for acc in missing_meta:
        print(f"  {acc}")

if unknown_admin2:
    print(f"\nWARNING: {len(unknown_admin2)} admin2 values not in aliases file:")
    for name in sorted(unknown_admin2):
        print(f"  {name}")

# Write output
out_path = f"{args.prefix}.aln"
writeFASTA(out_seqs, out_path)
print(f"\nWrote {out_path}")
