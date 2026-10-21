
from pathlib import Path
import argparse

def count_sequences(fasta_path):
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def filter_sequences(input_aln, output_aln, exclude_files):
    """Remove sequences whose ID or accession appears in any exclude file."""
    to_remove = set()
    for ef in exclude_files:
        if not Path(ef).exists():
            continue
        with open(ef) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                token = line.split()[0]
                if '|' in token:
                    to_remove.add(token.split('|')[0])
                    to_remove.add(token.split('|')[1].split('.')[0])
                else:
                    to_remove.add(token.split('.')[0])

    kept = 0
    removed = 0
    with open(output_aln, 'w') as out:
        with open(input_aln) as inp:
            writing = True
            for line in inp:
                if line.startswith('>'):
                    header = line[1:].strip()
                    fields = header.split('|')
                    seq_id = fields[0]
                    acc = fields[1].split('.')[0] if len(fields) > 1 else header.split('.')[0]
                    if seq_id in to_remove or acc in to_remove:
                        writing = False
                        removed += 1
                    else:
                        writing = True
                        kept += 1
                        out.write(line)
                elif writing:
                    out.write(line)
    print(f"Filtered: {kept} kept, {removed} removed")
    return kept

def main():
    parser = argparse.ArgumentParser(description="Filter sequences from an alignment file based on exclude files.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the alignment fasta file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output alignment fasta file.")
    parser.add_argument("-e", "--exclude-files", type=str, nargs='*', default=[], help="One or more text files with accession IDs to exclude (one per line, # comments allowed).")

    args = parser.parse_args()

    filter_sequences(args.input, args.output, args.exclude_files)

if __name__ == "__main__":
    main()

