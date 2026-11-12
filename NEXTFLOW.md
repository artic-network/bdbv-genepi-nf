# Running the pipeline with Nextflow

`main.nf` runs download → clean → QC → ADAR handling → temporal rooting → outlier filtering → final tree → health checks, with each step as its own process. Every process runs in the `bdbv-pipeline` Docker image by default.

## Setup

Install Nextflow (>= 23.04) and Docker. The pipeline pulls `articnetworkorg/bdbv-pipeline` from Docker Hub automatically; no build is needed.

To use a different image tag, e.g. your own build, set `--wf.container_version`:

```bash
docker build -t articnetworkorg/bdbv-pipeline:dev .
nextflow run . --wf.container_version dev ...
```

Other profiles: `-profile singularity`, or `-profile conda` to use `environment.yml` directly.

## Initial run

Downloads from Pathoplexus (unless `--fasta` and `--metadata` are given), builds the tree, auto-removing outliers with |z| ≥ 3:

```bash
nextflow run . --outdir data/YYYYMMDD --scope drc --min_date 2026-05-01 --adar_mode mask --outlier_removal_threshold 3
```

To use files you have already downloaded, add `--fasta data/YYYYMMDD/aligned.fasta --metadata data/YYYYMMDD/metadata.tsv`. Output file names start with `--prefix` (default `bdbv_<run_name>`, e.g. `bdbv_20261026`).

## Update run

As the initial run, plus `--previous_data` to run health checks against the previous dataset:

```bash
nextflow run . --outdir data/YYYYMMDD --previous_data data/YYYYMMDD_previous --outlier_removal_threshold 3
```

The report is written to `health_check_report.txt`. Unlike `bin/run.py`, the pipeline does not edit files in `data/`. It writes suggested additions, which you review and merge yourself:

- `exclude_pending.additions.txt` — new `duplicate_sample` flags to add to `data/exclude_pending.txt`
- `retractions.additions.txt` — retracted accessions to add to your `retractions.txt`

## Manual review point and resuming

There is a place where you may want to stop and review. Resume with `-resume` so completed steps are reused from cache, keeping the same `--outdir`, `--prefix` and input parameters.

**1. Outliers (after temporal rooting).** If you do not set `--outlier_removal_threshold`, the pipeline stops after temporal rooting. Review the RTT plot and `<prefix>.ADAR_masked.outliers.txt` in `temporal_root/`, copy the outliers file, edit the copy, then resume:

```bash
nextflow run . <same parameters> -resume --outlier_file my_outliers.txt
```

Health checks (update runs) are the final step. Review `health_check_report.txt` and update the exclusion lists. If something must be excluded, rerun with the updated `data/exclude_pending.txt` (the tree steps rerun because the `pre_filter` input changed).

## Parameters

| Parameter | Description |
|-----------|-------------|
| `--outdir` | Output folder (default `output`) |
| `--prefix` | Output file name prefix (default `bdbv_<run_name>`, e.g. `bdbv_20261026`) |
| `--run_name` | Label used in the default prefix and download file names (default: today, `YYYYMMDD`) |
| `--fasta`, `--metadata` | Use these instead of downloading |
| `--scope`, `--organism` | Download scope (`drc`, `drc,kenya`, `outbreak`) and LAPIS organism |
| `--max_date` | Latest collection date to download (`YYYY-MM-DD`; default no upper limit) |
| `--mafft` | Run MAFFT before cleaning |
| `--min_date` | Remove sequences dated before this (`YYYY-MM-DD`) |
| `--trim_start`, `--trim_end` | Alignment positions to keep, 1-based inclusive (default `47` and `18851`) |
| `--mask_sites`, `--mask_sites_file` | 1-based positions to mask to N (comma-separated list, or a file with one per line). The file defaults to `data/mask_sites.txt`; use `--mask_sites_file false` to mask none |
| `--output_intermediates` | Publish intermediate alignments (default `true`; `--output_intermediates false` to skip) |
| `--exclude_file`, `--exclude_pending` | Manual exclusions (default `data/exclude.txt`) and pending-review exclusions |
| `--model` | Skip ModelFinder and use this model |
| `--adar_mode` | `mask` (default), `exclude_tc` or `exclude_all` |
| `--seed` | IQ-TREE seed (default 42) |
| `--outlier_removal_threshold` | Auto-remove outliers with \|z\| ≥ N |
| `--outlier_file` | Reviewed outlier file (continues after review point 1) |
| `--previous_data` | Previous run folder, enables health checks |
| `--full_health_check` | Health checks on all genomes, not just new ones |

See `nextflow_schema.json` for all options.
