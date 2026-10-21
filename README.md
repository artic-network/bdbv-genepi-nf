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

### 2. Align with MAFFT

```bash
mafft --auto --thread -1 data/YYYYMMDD/aligned.fasta > data/YYYYMMDD/aligned_mafft.fasta
```

### 3. Run the pipeline

```bash
python scripts/run_pipeline.py \
    -f data/YYYYMMDD/aligned_mafft.fasta \
    -m data/YYYYMMDD/metadata.tsv \
    -p data/YYYYMMDD/prefix \
    --trim-start 1 --trim-end 18420 \
    --min-date 2026-05-01 \
    --exclude-file data/exclude_identity.txt \
    --exclude-quality data/exclude_quality.txt
```

The pipeline runs 12 steps:

1. **clean** — Rename headers, trim, deduplicate versions, apply date filter
2. **qc_basic** — Genome-level QC report
3. **pre_filter** — Remove excluded sequences
4. **iqtree1** — ML tree + ModelFinder (or `--model GTR+F+R3`) + ancestral reconstruction
5. **raccoon** — ADAR site identification
6. **adar_mask** — Mask ADAR-edited sites (default) or exclude affected genomes
7. **qc_phylo** — Phylogenetic QC report
8. **iqtree2** — Rebuild tree on masked alignment
9. **temporal_root** — Best-fitting temporal root + outlier detection

**← Stops here for manual outlier review →**

10. **filter_outliers** — Remove confirmed outliers
11. **iqtree3** — Final ML tree with aBayes branch support
12. **final_root** — Final temporal rooting

Resume after review:

```bash
python scripts/run_pipeline.py ... --skip-to filter_outliers
```

Or auto-remove outliers with z ≥ 3 and run straight through:

```bash
python scripts/run_pipeline.py ... --outlier-removal-threshold 3
```

### Key flags

| Flag | Description |
|------|-------------|
| `--model MODEL` | Skip ModelFinder, use this model directly |
| `--adar-mode` | `mask` (default), `exclude_tc`, or `exclude_all` |
| `--skip-to STEP` | Resume from a specific step |
| `--stop-after STEP` | Stop after a specific step |
| `--seed N` | IQ-TREE random seed (default: 42) |
