# BDBV Phylogenetic Pipeline

Automated pipeline for downloading, cleaning, and building maximum-likelihood phylogenies from Bundibugyo ebolavirus (BDBV) genome sequences hosted on [Pathoplexus](https://pathoplexus.org).

## Setup

### Dependencies

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate bdbv-pipeline
```

[IQ-TREE 2/3](http://www.iqtree.org/) must be installed separately and available on `PATH` as `iqtree3`.

### Data files

- **`data/aliases.csv`** — Maps observed admin2 location names to canonical names. This file needs to be updated periodically as new sequences with novel location spellings are submitted.
- **`data/exclude_identity.txt`** — Accessions to exclude from all analyses (duplicates, non-consortium, withdrawn).
- **`data/exclude_quality.txt`** — Accessions to exclude from phylogenetic analyses only (low quality).

## Usage

### 1. Download sequences

```bash
python scripts/download_pathoplexus.py --scope drc -o data/YYYYMMDD
```

`--scope` options:
- `drc` — DRC sequences only (applies Pathogen Genomic Lab group filter)
- `drc,kenya` — DRC + Kenya
- `outbreak` — All current outbreak countries (DRC, Kenya, Uganda)

### 2. Run the pipeline

Pathoplexus returns pre-aligned sequences, so MAFFT is not needed by default. Use `--mafft` if your input sequences are unaligned.

#### Automatic mode (recommended for routine updates)

Runs start to finish, auto-removing outliers with z-score ≥ 3:

```bash
python scripts/run_pipeline.py \
    -f data/YYYYMMDD/aligned.fasta \
    -m data/YYYYMMDD/metadata.tsv \
    -p data/YYYYMMDD/prefix \
    --min-date 2026-05-01 \
    --exclude-file data/exclude_identity.txt \
    --exclude-quality data/exclude_quality.txt \
    --adar-mode mask \
    --outlier-removal-threshold 3
```

#### Manual review mode

Omit `--outlier-removal-threshold` to stop after temporal rooting for manual outlier review:

```bash
python scripts/run_pipeline.py \
    -f data/YYYYMMDD/aligned.fasta \
    -m data/YYYYMMDD/metadata.tsv \
    -p data/YYYYMMDD/prefix \
    --min-date 2026-05-01 \
    --exclude-file data/exclude_identity.txt \
    --exclude-quality data/exclude_quality.txt \
    --adar-mode mask
```

The pipeline will stop after step 9 and print paths to the RTT plot and outlier file. Review and edit `outliers.txt`, then resume:

```bash
python scripts/run_pipeline.py ... --skip-to filter_outliers --adar-mode mask
```

### Pipeline steps

1. **mafft** *(optional, `--mafft`)* — MAFFT alignment
2. **clean** — Rename headers, trim, deduplicate versions, apply date filter
3. **qc_basic** — Genome-level QC report
4. **pre_filter** — Remove excluded sequences
5. **iqtree1** — ML tree + ModelFinder (or `--model MODEL`) + ancestral reconstruction
6. **raccoon** — ADAR site identification
7. **adar_mask** — Mask ADAR-edited sites or exclude affected genomes
8. **qc_phylo** — Phylogenetic QC report
9. **iqtree2** — Rebuild tree on masked alignment
10. **temporal_root** — Best-fitting temporal root + outlier detection
11. **filter_outliers** — Remove outliers (auto or manual)
12. **iqtree3** — Final ML tree with aBayes branch support
13. **final_root** — Final temporal rooting

### Key flags

| Flag | Description |
|------|-------------|
| `--outlier-removal-threshold N` | Auto-remove outliers with \|z\| ≥ N and run straight through (e.g. `3`) |
| `--model MODEL` | Skip ModelFinder, use this model directly (e.g. `GTR+F+R3`) |
| `--mafft` | Run MAFFT alignment before cleaning |
| `--adar-mode` | `mask` (default), `exclude_tc`, or `exclude_all` |
| `--skip-to STEP` | Resume from a specific step |
| `--stop-after STEP` | Stop after a specific step |
| `--seed N` | IQ-TREE random seed (default: 42) |
| `--min-date YYYY-MM-DD` | Remove sequences before this date; same-year sequences are flagged |
