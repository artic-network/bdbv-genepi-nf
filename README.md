# BDBV Phylogenetic Pipeline

Automated pipeline for downloading, cleaning, and building maximum-likelihood phylogenies from Bundibugyo ebolavirus (BDBV) genome sequences hosted on [Pathoplexus](https://pathoplexus.org).

## Setup

### Dependencies

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate bdbv-pipeline
```

`environment.yml` includes IQ-TREE 3, which conda can't install on every platform. If it fails on yours, remove that line and install [IQ-TREE 2/3](http://www.iqtree.org/) separately, available on `PATH` as `iqtree3`. The Docker image used by the Nextflow pipeline always includes it.

### Data files

- **`data/aliases.csv`** — Maps observed admin2 location names to canonical names. This file needs to be updated periodically as new sequences with novel location spellings are submitted.
- **`data/exclude_identity.txt`** — Accessions to exclude from all analyses (duplicates, non-consortium, withdrawn).
- **`data/exclude_quality.txt`** — Accessions to exclude from phylogenetic analyses only (low quality).
- **`data/exclude_pending.txt`** — Accessions temporarily excluded while under review (e.g. flagged `duplicate_sample` by health checks).
- **`data/exclude_pending.txt`** — Accessions temporarily excluded while under review.

## Running with Nextflow

The Nextflow pipeline (`main.nf`) runs download → tree building → health checks → BEAST XMLs, with each step as a separate process. By default every step runs in Docker.

### Setup

Install [Nextflow](https://www.nextflow.io/) (≥ 23.04) and Docker, then build the image from this folder:

```bash
docker build -t bdbv-pipeline:latest .
```

On Apple Silicon add `--platform linux/amd64`. To use conda instead of Docker, add `-profile conda` to the commands below.

### Initial run

Downloads the latest sequences from Pathoplexus and builds the trees and BEAST XMLs:

```bash
nextflow run . --outdir data/YYYYMMDD --min_date 2026-05-01 --outlier_removal_threshold 3
```

To use files you've already downloaded inst

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

### 3. Generate BEAST XMLs

After the pipeline completes, generate BEAST XML files from the final alignment using template files:

```bash
python scripts/generate_beast_xml.py \
    -a data/YYYYMMDD/prefix.ADAR_masked.final.aln \
    -o data/YYYYMMDD \
    -p prefix.n978 \
    --template-egc templates/GTR_EGC.template.xml \
    --template-sg templates/GTR_SG.template.xml \
    --template-lgc templates/GTR_LGC.template.xml \
    --anchor-date 2025-11-08
```

This swaps the alignment and taxa into the template, updates SkyGrid parameters (weekly grid intervals from anchor date to the latest sample), adjusts operator weights for the dataset size, and sets MCMC chain length / logging intervals.

| Flag | Description |
|------|-------------|
| `--template-egc` | Exponential growth coalescent template |
| `--template-sg` | SkyGrid coalescent template (requires `--anchor-date`) |
| `--template-lgc` | Logistic growth coalescent template |
| `--anchor-date` | SkyGrid anchor date — grid extends from this date to the latest sample. We currently use `2025-11-08`. |
| `--chain-length` | MCMC chain length (default: 500M) |
| `--log-every` | Log sampling interval (default: 50,000) |

Templates are in `templates/` with placeholder taxa and N-filled sequences. The script replaces these with the real data.

### 4. One-command update

`run_update.py` chains download → pipeline → health checks → BEAST XMLs into a single command:

```bash
python scripts/run_update.py \
    --previous-data data/20260915 \
    --template-sg templates/GTR_SG.template.xml \
    --template-egc templates/GTR_EGC.template.xml \
    --anchor-date 2025-11-08
```

This automatically:
- Downloads the latest sequences from Pathoplexus
- Checks for new genomes (exits early if count is unchanged)
- Merges exclude lists (identity + quality + pending)
- Runs the full pipeline with auto-outlier removal
- Runs health checks comparing to the previous dataset
- Auto-appends `duplicate_sample` flags to `data/exclude_pending.txt`
- Generates BEAST XMLs

Use `--skip-download`, `--skip-pipeline`, or `--skip-beast-xmls` to skip individual steps.

### 5. Health checks

`health_checks.py` compares current and previous pipeline runs:

```bash
python scripts/health_checks.py --current data/20261012 --previous data/20260915 --prefix Ituri2026.DRC_trimmed
```

Per-sequence flags: `high_adar`, `high_Ns`, `low_coverage`, `near_outlier`, `many_gaps`, `long_branch`, `duplicate_sample`, `phylo_conflict`.

Dataset-level checks: outlier count, clock rate shift, R², ADAR branch count, genome count, version changes/retractions.

## Environment notes

- **raccoon** requires `baltic<1` (pinned in `environment.yml`)
- **IQ-TREE** is in `environment.yml` and so in the Docker image, but is not available via conda for all platforms. If it is missing when running the scripts directly, install it separately
