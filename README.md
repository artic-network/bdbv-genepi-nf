# BDBV Phylogenetic Pipeline

Automated pipeline for downloading, cleaning, and building maximum-likelihood phylogenies from Bundibugyo ebolavirus (BDBV) genome sequences hosted on [Pathoplexus](https://pathoplexus.org).

## Setup

### Dependencies

Create the conda environment:

```bash
mamba env create -f environment.yml
```

All pipeline commands should be run with `mamba run -n bdbv` (or `conda run -n bdbv`):

```bash
mamba run -n bdbv python bin/run.py ...
```

[IQ-TREE 2/3](http://www.iqtree.org/) must be installed separately and available on `PATH` as `iqtree3`.

### Data files

- **`data/aliases.csv`** — Maps observed admin2 location names to canonical health zone names. Update as new location spellings appear.
- **`data/exclude.txt`** — Manual exclusions (accessions to always exclude from analysis).
- **`data/exclusion_log.tsv`** — Audit trail documenting exclude list changes.
- **`data/health_area_zone_province_slim.csv`** — Health area → zone → province lookup for location validation.

Per-run auto-generated files (not checked into the repo):
- **`exclude_auto.txt`** — Resolved duplicates, unresolved duplicates (temporarily excluded), and quality failures.
- **`exclude_merged.txt`** — Combined manual + auto exclusions, used by the pipeline.

## Usage

### `run.py` — single entry point

Minimal (DRC consortium only, no comparison):

```bash
mamba run -n bdbv python bin/run.py \
    --outdir data/20261104_drc \
    --prefix BDBV2026.drc \
    --scope drc \
    --model GTR+F+R3 \
    --trim-start 47 --trim-end 18851 \
    --mask-sites-file data/mask_sites.txt \
    --seed 42
```

Multi-country outbreak scope:

```bash
mamba run -n bdbv python bin/run.py \
    --outdir data/20261104_outbreak \
    --prefix BDBV2026.outbreak \
    --scope outbreak \
    --model GTR+F+R3 \
    --trim-start 47 --trim-end 18851 \
    --mask-sites-file data/mask_sites.txt \
    --seed 42
```

Full run with linelist, previous data, and BEAST XMLs:

```bash
mamba run -n bdbv python bin/run.py \
    --previous-data data/20261027_drc \
    --outdir data/20261104_drc \
    --prefix BDBV2026.drc \
    --scope drc \
    --model GTR+F+R3 \
    --trim-start 47 --trim-end 18851 \
    --mask-sites-file data/mask_sites.txt \
    --seed 42 \
    --linelist data/CONFIDENTIAL-linelist.xlsx \
    --aliases-file data/aliases.csv \
    --hz-lookup data/health_area_zone_province_slim.csv \
    --template-sg templates/GTR_SG_MAT.template.xml
```

**Steps** (each logged with timestamps to `run_log.txt`):

1. **Download** — Fetch latest sequences from Pathoplexus
2. **Count check** — Compare to previous run; exit early if no new genomes *(requires `--previous-data`)*
3. **Metadata quality report** — Linelist linking, duplicate detection, location/date validation, sequence quality checks
4. **Merge exclude lists** — Combine `exclude.txt` + `exclude_auto.txt` → `exclude_merged.txt`
5. **Pipeline** — Clean, ADAR masking, tree building, temporal rooting, outlier removal
6. **Health checks** — Per-sequence flags and dataset-level checks *(requires `--previous-data`)*
7. **BEAST XMLs** — Generate XML files from templates *(requires `--template-sg` or `--template-egc`)*
8. **Run comparison** — Side-by-side diff with previous run *(requires `--previous-data`)*
9. **Summary** — List output files

| Flag | Description |
|------|-------------|
| `--previous-data` | Previous pipeline output directory (enables comparison and health checks) |
| `--outdir` | Output directory (default: `data/YYYYMMDD`) |
| `--prefix` | Output file prefix |
| `--scope` | `drc` (consortium only), `outbreak` (all groups, all outbreak countries), or comma-separated country list |
| `--model` | IQ-TREE substitution model (e.g. `GTR+F+R3`; default: ModelFinder) |
| `--trim-start` | First alignment position to keep, 1-based inclusive (e.g. `47`) |
| `--trim-end` | Last alignment position to keep, 1-based inclusive (e.g. `18851`) |
| `--mask-sites` | Comma-separated 1-based positions to mask to N (e.g. `5246,5249,5251,5252,5255`) |
| `--mask-sites-file` | File with 1-based positions to mask, one per line (`#` comments allowed) |
| `--outlier-threshold` | Auto-remove outliers with \|z\| ≥ this value (default: `3.0`) |
| `--seed` | IQ-TREE random seed (default: `42`) |
| `--linelist` | INRB linelist Excel file (enables duplicate resolution and mismatch checking) |
| `--exclude` | Manual exclude list (default: `data/exclude.txt`) |
| `--min-date` | Earliest collection date to download (default: `2025-01-01`) |
| `--adar-mode` | `mask` (default), `exclude_tc`, or `exclude_all` |
| `--skip-download` | Use existing files in `--outdir` |
| `--skip-pipeline` | Skip pipeline step |
| `--skip-beast-xmls` | Skip BEAST XML generation |
| `--full-check` | Run health checks on all genomes, not just new ones |
| `--template-sg` | SkyGrid template XML |
| `--template-egc` | EGC template XML |
| `--anchor-date` | SkyGrid anchor date (default: `2025-11-08`) |
| `--slack-webhook` | Slack webhook URL for notifications |

### Inputs and outputs

**Inputs:**
- Pathoplexus API (sequences + metadata, downloaded automatically)
- INRB linelist Excel (for duplicate resolution and metadata validation)
- `data/exclude.txt` (manual exclusions)
- `data/aliases.csv` (location name mappings)
- `data/health_area_zone_province_slim.csv` (location hierarchy)
- BEAST template XMLs (optional, for XML generation)

**Outputs (in `--outdir`):**

| File | Description |
|------|-------------|
| `aligned.fasta` | Raw downloaded alignment |
| `metadata.tsv` | Raw downloaded metadata |
| `metadata_quality_report.txt` | Human-readable quality report for the team |
| `metadata_quality_report.csv` | Machine-readable quality issues (one row per issue) |
| `non_consortium_sequences.csv` | Non-consortium sequences separated out |
| `exclude_auto.txt` | Auto-generated exclusions (duplicates + quality) |
| `exclude_merged.txt` | Combined manual + auto exclusions |
| `{prefix}.aln` | Cleaned alignment |
| `{prefix}.ADAR_masked.aln` | ADAR-masked alignment |
| `{prefix}.ADAR_masked.adar_report.txt` | ADAR editing details |
| `{prefix}.ADAR_masked.final.aln` | Final alignment (after outlier removal) |
| `{prefix}.ADAR_masked.final.treefile` | Final ML tree |
| `{prefix}.ADAR_masked.final.rtt_summary.txt` | Clock rate and R² |
| `{prefix}.ADAR_masked.final.outliers.txt` | Removed outliers with z-scores |
| `health_check_report.txt` | Per-sequence flags and dataset-level checks |
| `run_comparison.txt` | Diff vs previous run |
| `run_log.txt` | Timestamped log of all steps |
| `*.xml` | BEAST XML files (if templates provided) |

## Metadata quality report

`metadata_quality_report.py` validates consortium sequences against the INRB linelist and checks sequence quality.

**Linelist linking** tries three keys in order:
1. `specimenCollectorSampleId` → linelist `Sample ID`
2. `submissionId` → linelist `Sample ID`
3. `accessionVersion` → linelist `Pathoplexus ID`

**Report sections:**

*Part 1 — Needs input:*
1. Not in linelist (has an ID but no linelist match)
2. Unresolved duplicates (neither side in linelist; temporarily keeping higher-completeness one)
3. Missing `specimenCollectorSampleId`
4. Linelist mismatches (date or location differs between Pathoplexus and linelist)
5. Location issues (health area instead of zone, province instead of zone, missing)
6. Date issues (incomplete or missing collection dates)

*Auto-excluded:*
- Completeness below 70%
- Non-ADAR mismatches above 10%

*Part 2 — Already resolved:*
- Resolved duplicates with KEEP/REVOKE and resolution method (linelist match)

**Duplicate resolution order:**
1. Linelist match (Tier 1: explicit Pathoplexus ID match; Tier 2: one side linked, the other not)
2. If neither side is in the linelist → unresolved (temporarily keep higher completeness, flag for team)

Non-consortium sequences (any group other than `Pathogen Genomic Lab`) are separated into `non_consortium_sequences.csv` and excluded from the team-facing report.

## Individual scripts

### Download

```bash
mamba run -n bdbv python bin/download_pathoplexus.py --scope outbreak -o data/YYYYMMDD
```

`--scope` options:
- `drc` — DRC consortium sequences only (`Pathogen Genomic Lab` group filter)
- `outbreak` — All groups, all outbreak countries (DRC, Kenya, Nigeria, Uganda)

Defaults to `--min-date 2025-01-01` to exclude historical sequences.

### Pipeline (advanced)

`run_pipeline.py` is the inner engine called by `run.py`. Use it directly only when you need fine-grained step control:

```bash
mamba run -n bdbv python bin/run_pipeline.py \
    -f data/YYYYMMDD/aligned.fasta \
    -m data/YYYYMMDD/metadata.tsv \
    -p prefix \
    -d data/YYYYMMDD \
    --exclude-file data/YYYYMMDD/exclude_merged.txt \
    --adar-mode mask \
    --trim-start 47 --trim-end 18851 \
    --seed 42 \
    --outlier-removal-threshold 3
```

Omit `--outlier-removal-threshold` for manual outlier review. The pipeline stops after temporal rooting; review `outliers.txt`, then resume:

```bash
mamba run -n bdbv python bin/run_pipeline.py \
    -f data/YYYYMMDD/aligned.fasta \
    -m data/YYYYMMDD/metadata.tsv \
    -p prefix \
    -d data/YYYYMMDD \
    --adar-mode mask \
    --trim-start 47 --trim-end 18851 \
    --seed 42 \
    --outlier-removal-threshold 3 \
    --skip-to filter_outliers
```

#### Pipeline steps

1. **mafft** *(optional)* — MAFFT alignment
2. **clean** — Rename headers, trim, mask sites, deduplicate versions
3. **qc_basic** — Genome-level QC report
4. **pre_filter** — Remove excluded sequences
5. **iqtree1** — ML tree + ModelFinder + ancestral reconstruction
6. **raccoon** — ADAR site identification
7. **adar_mask** — Mask ADAR-edited sites or exclude affected genomes
8. **qc_phylo** — Phylogenetic QC report
9. **iqtree2** — Rebuild tree on masked alignment
10. **temporal_root** — Best-fitting temporal root + outlier detection
11. **filter_outliers** — Remove outliers (auto or manual)
12. **iqtree3** — Final ML tree with aBayes branch support
13. **final_root** — Final temporal rooting

| Flag | Description |
|------|-------------|
| `--outlier-removal-threshold N` | Auto-remove outliers with \|z\| ≥ N (e.g. `3`) |
| `--trim-start N` | First alignment position to keep, 1-based inclusive (e.g. `47`) |
| `--trim-end N` | Last alignment position to keep, 1-based inclusive (e.g. `18851`) |
| `--mask-sites` | Comma-separated 1-based positions to mask to N (e.g. `5246,5249,5251,5252,5255`) |
| `--mask-sites-file` | File with 1-based positions to mask, one per line (`#` comments allowed) |
| `--adar-mode` | `mask` (default), `exclude_tc`, or `exclude_all` |
| `--skip-to STEP` | Resume from a specific step |
| `--stop-after STEP` | Stop after a specific step |
| `--seed N` | IQ-TREE random seed (default: 42) |

**Mask sites file format** (one 1-based position per line):

```
# Primer-trimming artifact at VP40 amplicon 13/14 junction
5246
5249
5251
5252
5255
```

### BEAST XMLs

```bash
mamba run -n bdbv python bin/generate_beast_xml.py \
    -a data/YYYYMMDD/prefix.ADAR_masked.final.aln \
    -o data/YYYYMMDD \
    -p prefix.nNNN \
    --template-sg templates/GTR_SG.template.xml \
    --template-egc templates/GTR_EGC.template.xml \
    --anchor-date 2025-11-08
```

Templates are in `templates/` with placeholder taxa. The script swaps in the real alignment, updates SkyGrid parameters (weekly grid intervals from anchor date to latest sample), and adjusts operator weights for dataset size.

Two SkyGrid templates are available:
- **`GTR_SG.template.xml`** — Standard SkyGrid. Works with any BEAST X version.
- **`GTR_SG_MAT.template.xml`** — SkyGrid using mutation-annotated trees. Requires a MAT-specific release of BEAST X.

### Health checks

```bash
mamba run -n bdbv python bin/health_checks.py --current data/20261029 --previous data/20261012 --prefix Ituri2026.DRC_trimmed
```

Per-sequence flags (new genomes by default, `--full` for all): `high_adar`, `high_Ns`, `low_coverage`, `near_outlier`, `many_gaps`, `long_branch`, `duplicate_sample`, `phylo_conflict`.

Dataset-level checks: outlier count, clock rate shift, R², ADAR branch count, genome count, version changes/retractions.

### Run comparison

```bash
mamba run -n bdbv python bin/compare_runs.py --current data/20261029 --previous data/20261012 --prefix Ituri2026.DRC_trimmed
```

Side-by-side diff of two pipeline runs: sequence counts, added/removed accessions, exclude list changes, clock rate, ADAR summary, metadata differences.

## Environment notes

- **raccoon** requires `baltic<1` (pinned in `environment.yml`)
- **IQ-TREE** must be installed separately (not available via conda for all platforms)
