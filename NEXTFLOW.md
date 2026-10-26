# Running the pipeline with Nextflow

`main.nf` runs download → clean → QC → ADAR handling → temporal rooting → outlier filtering → final tree → health checks → BEAST XMLs, with each step as its own process. Every process runs in the `bdbv-pipeline` Docker image by default.

## Setup

Install Nextflow (>= 23.04) and Docker. The pipeline pulls `articnetworkorg/bdbv-pipeline` from Docker Hub automatically; no build is needed.

To use a different image tag, e.g. your own build, set `--wf.container_version`:

```bash
docker build -t articnetworkorg/bdbv-pipeline:dev .
nextflow run . --wf.container_version dev ...
```

Other profiles: `-profile singularity`, or `-profile conda` to use `environment.yml` directly.

## Initial run

Downloads from Pathoplexus (unless `--fasta` and `--metadata` are given), builds the tree and generates BEAST XMLs, auto-removing outliers with |z| ≥ 3:

```bash
nextflow run . --outdir data/YYYYMMDD --scope drc --min_date 2026-05-01 --adar_mode mask --outlier_removal_threshold 3
```

To use files you have already downloaded, add `--fasta data/YYYYMMDD/aligned.fasta --metadata data/YYYYMMDD/metadata.tsv`. Output file names start with `--prefix` (default `bdbv-pipeline`).

## Update run

As the initial run, plus `--previous_data` to run health checks against the previous dataset:

```bash
nextflow run . --outdir data/YYYYMMDD --previous_data data/YYYYMMDD_previous --outlier_removal_threshold 3
```

The report is written to `health_check_report.txt`. Unlike `run_update.py`, the pipeline does not edit files in `data/`. It writes suggested additions, which you review and merge yourself:

- `exclude_pending.additions.txt` — new `duplicate_sample` flags to add to `data/exclude_pending.txt`
- `retractions.additions.txt` — retracted accessions to add to your `retractions.txt`

## Manual review points and resuming

There are two places where you may want to stop and review. Resume with `-resume` so completed steps are reused from cache, keeping the same `--outdir`, `--prefix` and input parameters.

**1. Outliers (after temporal rooting).** If you do not set `--outlier_removal_threshold`, the pipeline stops after temporal rooting. Review the RTT plot and `<prefix>.ADAR_masked.outliers.txt` in the output folder, copy the outliers file, edit the copy, then resume:

```bash
nextflow run . <same parameters> -resume --outlier_file my_outliers.txt
```

**2. Health checks (before BEAST XMLs).** Add `--stop_after_health_checks` (with `--previous_data`) to stop after health checks. Review the report and update the exclusion lists. If something must be excluded, rerun with the updated `data/exclude_pending.txt` (the tree steps rerun because the `pre_filter` input changed). If nothing needs excluding, resume without the flag to generate the XMLs:

```bash
nextflow run . <same parameters> -resume
```

## Parameters

| Parameter | Description |
|-----------|-------------|
| `--outdir` | Output folder (default `output`) |
| `--prefix` | Output file name prefix (default `bdbv-pipeline`) |
| `--fasta`, `--metadata` | Use these instead of downloading |
| `--scope`, `--organism` | Download scope (`drc`, `drc,kenya`, `outbreak`) and LAPIS organism |
| `--mafft` | Run MAFFT before cleaning |
| `--min_date`, `--trim_start`, `--trim_end` | Cleaning options |
| `--model` | Skip ModelFinder and use this model |
| `--adar_mode` | `mask` (default), `exclude_tc` or `exclude_all` |
| `--seed` | IQ-TREE seed (default 42) |
| `--outlier_removal_threshold` | Auto-remove outliers with \|z\| ≥ N |
| `--outlier_file` | Reviewed outlier file (continues after review point 1) |
| `--previous_data` | Previous run folder, enables health checks |
| `--full_health_check` | Health checks on all genomes, not just new ones |
| `--stop_after_health_checks` | Stop for review point 2 |
| `--beast_xml false` | Skip BEAST XML generation |
| `--template_egc`, `--template_sg`, `--template_lgc` | BEAST templates (EGC and SG by default, LGC off) |
| `--anchor_date` | SkyGrid anchor date (default `2025-11-08`) |

See `nextflow_schema.json` for all options.
