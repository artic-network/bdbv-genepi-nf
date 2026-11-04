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
- **`data/mask_sites.txt`** — Alignment positions to mask to N before phylogenetic analysis. Currently contains 5 sites (5246, 5249, 5251, 5252, 5255) in the VP40 gene at the amplicon 13/14 junction. These are primer-trimming artifacts from an Oct 2026 sequencing batch that was assembled with an incorrect BED file, causing untrimmed primer bases to be incorporated into consensus sequences as systematic miscalls. Site 5246 is the most homoplasic site in the genome (60 inferred mutation events on the tree, 3.5× more than any other site). See `docs/site_5246_primer_artifact.md` in the main project repo for the full diagnostic report.

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

Most users should use `run.py` (see above), which calls these scripts in order. Use individual scripts when you need to rerun a single step, debug, or work outside the standard workflow.

All commands should be prefixed with `mamba run --live-stream -n bdbv`.

### `download_pathoplexus.py` — Download sequences

Fetches aligned sequences and metadata from the Pathoplexus LAPIS API.

```bash
python bin/download_pathoplexus.py --scope drc -o /path/to/output/dir
```

| Flag | Description |
|------|-------------|
| `--scope` | `drc` (consortium only, filters to `Pathogen Genomic Lab` group), `outbreak` (all groups, all outbreak countries: DRC, Kenya, Nigeria, Uganda), or comma-separated country list |
| `-o` | Output directory |
| `--min-date` | Earliest collection date (default: `2025-01-01`) |
| `--organism` | Pathoplexus organism slug (default: `ebola-bdbv`) |

**When to use standalone:** Testing download without running the full pipeline, or downloading for a custom scope.

### `metadata_quality_report.py` — Metadata quality report

Validates sequences against the INRB linelist. Detects duplicates, linelist mismatches, location issues, date issues, and low-quality sequences. Outputs an actionable report for the team and auto-generates an exclude list.

```bash
python bin/metadata_quality_report.py \
    -m /path/to/metadata.tsv \
    -f /path/to/aligned.fasta \
    -a data/aliases.csv \
    -l /path/to/linelist.xlsx \
    -s data/health_area_zone_province_slim.csv \
    -o /path/to/output/metadata_quality_report.csv \
    --report-txt /path/to/output/metadata_quality_report.txt \
    --exclude-out /path/to/output/exclude_auto.txt \
    --non-consortium-out /path/to/output/non_consortium_sequences.csv
```

| Flag | Description |
|------|-------------|
| `-m` | Pathoplexus metadata TSV (required) |
| `-f` | FASTA alignment (required, for completeness/quality checks) |
| `-a` | Aliases CSV (required) |
| `-l` | INRB linelist Excel (required for duplicate resolution; omit to skip linelist checks) |
| `-s` | Health area/zone/province CSV (for location validation) |
| `-o` | Output CSV (one row per issue) |
| `--report-txt` | Human-readable report |
| `--exclude-out` | Auto-generated exclude list (duplicates + quality failures) |
| `--non-consortium-out` | Non-consortium sequences CSV |

**When to use standalone:** Regenerating the quality report after linelist updates without rerunning the pipeline.

### `clean_up_pathoplexus.py` — Clean alignment

Standardizes FASTA headers, masks terminal positions and specific sites to N, deduplicates sequence versions, applies exclude lists, and filters by date. This is the first pipeline step that modifies the alignment.

```bash
python bin/clean_up_pathoplexus.py \
    -f /path/to/aligned.fasta \
    -m /path/to/metadata.tsv \
    -a data/aliases.csv \
    -p /path/to/output/prefix \
    --trim_start 47 --trim_end 18851 \
    --mask_sites_file data/mask_sites.txt \
    --min-date 2025-01-01
```

| Flag | Description |
|------|-------------|
| `-f` | Input FASTA (required) |
| `-m` | Pathoplexus metadata TSV (required) |
| `-p` | Output prefix (required) — produces `{prefix}.aln` |
| `-a` | Aliases CSV (default: `data/aliases.csv`) |
| `--trim_start` | First position to keep (1-based); earlier positions masked to N |
| `--trim_end` | Last position to keep (1-based); later positions masked to N |
| `--mask_sites` | Comma-separated 1-based positions to mask to N |
| `--mask_sites_file` | File with positions to mask (one per line, `#` comments allowed) |
| `-e` | One or more exclude files (accession IDs to remove) |
| `--min-date` | Remove sequences before this date (YYYY-MM-DD) |

**When to use standalone:** Testing trimming/masking parameters, or cleaning an alignment for use outside the pipeline.

### `genome_qc_report.py` — Basic genome QC

Reports per-sequence quality metrics (completeness, N fraction, gaps) on the cleaned alignment, before tree building.

```bash
python bin/genome_qc_report.py \
    -a /path/to/cleaned.aln \
    -m /path/to/metadata.tsv \
    -o /path/to/output/qc_report_basic.txt
```

**When to use standalone:** Quick quality check on an alignment without building a tree.

### `mask_adar_sites.py` — ADAR masking

Uses raccoon output to handle ADAR-edited sites. Three modes: `mask` (replace sites with N), `exclude_tc` (remove genomes with T>C clusters), `exclude_all` (remove genomes with any ADAR-consistent clusters).

```bash
python bin/mask_adar_sites.py \
    -a /path/to/alignment.aln \
    -r /path/to/phylo_flags.csv \
    -t /path/to/treefile \
    -s /path/to/statefile \
    --mode mask \
    -p /path/to/output/prefix
```

| Flag | Description |
|------|-------------|
| `-a` | Input alignment (required) |
| `-r` | Raccoon phylo_flags CSV (required) |
| `-t` | Newick tree file (required) |
| `-s` | IQ-TREE `.state` file (required) |
| `--mode` | `mask`, `exclude_tc`, or `exclude_all` (required) |
| `-p` | Output prefix — produces `{prefix}.aln` and `{prefix}.adar_report.txt` |

In `mask` mode: GP gene sites are masked to N across all genomes (whole-column); non-GP sites are masked per-branch, keeping the lowest T>C site per branch.

**When to use standalone:** Rerunning ADAR handling with a different mode without rebuilding the tree.

### `phylo_qc_report.py` — Phylogenetic QC

Reports per-sequence phylogenetic quality metrics using raccoon output: ADAR editing levels, branch-level SNP patterns, and conflicts between genome QC and tree placement.

```bash
python bin/phylo_qc_report.py \
    -f /path/to/phylo_flags.csv \
    -b /path/to/branch_snps.reconstruction.csv \
    -m /path/to/metadata.tsv \
    --basic-qc /path/to/qc_report_basic.txt \
    -o /path/to/output/qc_report_phylo.txt
```

**When to use standalone:** Investigating specific sequences flagged by raccoon.

### `temporal_root.py` — Temporal rooting and outlier detection

Finds the best-fitting temporal root (maximizes R² of root-to-tip regression), identifies outliers by residual z-score, and optionally iterates to refine the root after removing outliers.

```bash
python bin/temporal_root.py \
    -t /path/to/treefile \
    --iterate \
    --plot \
    -o /path/to/output/ \
    --prefix my_analysis
```

| Flag | Description |
|------|-------------|
| `-t` | Newick tree file (required) |
| `--iterate` | Iteratively re-root after removing outliers until stable |
| `--plot` | Generate RTT plot (PDF) |
| `-o` | Output directory |
| `--prefix` | Output file prefix |
| `-z` | Z-score threshold for outlier detection (default: 2.0) |
| `-d` | Tip label delimiter (default: `\|`) |
| `-f` | Date field index in tip labels (default: `-1`, last field) |

Outputs: `{prefix}.bestFitRoot.tree`, `{prefix}.outliers.txt`, `{prefix}.rtt_plot.pdf`, `{prefix}.rtt_summary.txt`.

**When to use standalone:** Rooting a tree built outside the pipeline, or testing different outlier thresholds.

### `run_pipeline.py` — Inner pipeline engine

Runs the 13-step ML tree pipeline. Called by `run.py` but can be used directly for fine-grained step control.

```bash
python bin/run_pipeline.py \
    -f /path/to/aligned.fasta \
    -m /path/to/metadata.tsv \
    -p prefix \
    -d /path/to/output/ \
    --exclude-file /path/to/exclude_merged.txt \
    --model GTR+F+R3 \
    --adar-mode mask \
    --trim-start 47 --trim-end 18851 \
    --mask-sites-file data/mask_sites.txt \
    --seed 42 \
    --outlier-removal-threshold 3
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
| `--model` | Substitution model (e.g. `GTR+F+R3`; default: ModelFinder) |
| `--adar-mode` | `mask` (default), `exclude_tc`, or `exclude_all` |
| `--skip-to STEP` | Resume from a specific step |
| `--stop-after STEP` | Stop after a specific step |
| `--seed N` | IQ-TREE random seed (default: 42) |

Omit `--outlier-removal-threshold` for manual outlier review. The pipeline stops after temporal rooting; review `outliers.txt`, then resume with `--skip-to filter_outliers`.

**Mask sites file format** (one 1-based position per line):

```
# Primer-trimming artifact at VP40 amplicon 13/14 junction
5246
5249
5251
5252
5255
```

**When to use standalone:** Resuming from a specific step (`--skip-to`), stopping early (`--stop-after`), or running with manual outlier review (omit `--outlier-removal-threshold`).

### `generate_beast_xml.py` — BEAST XML generation

Generates BEAST XML files from templates. Swaps in the real alignment, updates SkyGrid parameters (weekly grid intervals from anchor date to latest sample), and adjusts operator weights for dataset size.

```bash
python bin/generate_beast_xml.py \
    -a /path/to/prefix.ADAR_masked.final.aln \
    -o /path/to/output/ \
    -p prefix.nNNN \
    --template-sg templates/GTR_SG_MAT.template.xml \
    --anchor-date 2025-11-08
```

| Flag | Description |
|------|-------------|
| `-a` | Final alignment (required) |
| `-o` | Output directory |
| `-p` | Output prefix |
| `--template-sg` | SkyGrid template XML |
| `--template-egc` | EGC template XML |
| `--anchor-date` | SkyGrid anchor date (default: `2025-11-08`) |
| `--chain-length` | MCMC chain length (default: 100,000,000) |
| `--log-every` | Log sampling interval (default: 10,000) |

Two SkyGrid templates are available:
- **`GTR_SG.template.xml`** — Standard SkyGrid. Works with any BEAST X version.
- **`GTR_SG_MAT.template.xml`** — SkyGrid using mutation-annotated trees. Requires a MAT-specific release of BEAST X.

**When to use standalone:** Generating XMLs for a different alignment or with different chain settings.

### `health_checks.py` — Health checks

Runs per-sequence flags and dataset-level checks. Most useful when comparing two runs to detect regressions.

```bash
# With comparison to previous run:
python bin/health_checks.py \
    --current /path/to/20261104_drc \
    --previous /path/to/20261103_drc \
    --prefix BDBV2026.drc \
    -o /path/to/20261104_drc/health_check_report.txt

# Standalone (no comparison, dataset-level checks only):
python bin/health_checks.py \
    --current /path/to/20261104_drc \
    --prefix BDBV2026.drc \
    --full
```

| Flag | Description |
|------|-------------|
| `--current` | Current pipeline output directory (required) |
| `--previous` | Previous pipeline output directory (enables comparison) |
| `--prefix` | File prefix pattern |
| `--full` | Flag all genomes, not just new ones |
| `-o` | Save report to file (otherwise prints to stdout) |

**Per-sequence flags:** `high_adar`, `high_Ns`, `low_coverage`, `near_outlier`, `many_gaps`, `long_branch`, `duplicate_sample`, `phylo_conflict`.

**Dataset-level checks:** outlier count, clock rate shift, R², ADAR branch count, genome count, version changes/retractions.

Without `--previous`, comparison checks (rate shift, new/removed accessions, version changes, watchlist deltas) are skipped. The standalone checks are still useful as a sanity check, but most of the information is already in the pipeline output files (`rtt_summary.txt`, `outliers.txt`, `adar_report.txt`, `qc_report_*.txt`).

**When to use standalone:** Post-hoc comparison of two runs, or quick sanity check on a completed run without rerunning the pipeline.

### `compare_runs.py` — Run comparison

Side-by-side diff of two completed pipeline runs.

```bash
python bin/compare_runs.py \
    --current /path/to/20261104_drc \
    --previous /path/to/20261103_drc \
    --prefix BDBV2026.drc \
    --adar-mode mask \
    -o /path/to/20261104_drc/run_comparison.txt
```

Reports: sequence counts, added/removed accessions, exclude list changes, clock rate comparison, ADAR summary, metadata differences.

**When to use standalone:** Comparing any two runs after the fact, without rerunning the pipeline.

## Pipeline DAG (for Nextflow conversion)

Each row shows a pipeline step, the script it runs, its required inputs, its outputs, and external tool dependencies. File patterns use `{p}` for prefix and `{d}` for output directory.

| Step | Script | Inputs | Outputs | External tools |
|------|--------|--------|---------|----------------|
| download | `download_pathoplexus.py` | Pathoplexus API | `{d}/ebola-bdbv_aligned-nuc_*.fasta`, `{d}/ebola-bdbv_metadata_*.tsv` | — |
| metadata_qc | `metadata_quality_report.py` | `metadata.tsv`, `aligned.fasta`, `aliases.csv`, `linelist.xlsx` (optional) | `{d}/metadata_quality_report.csv`, `{d}/metadata_quality_report.txt`, `{d}/exclude_auto.txt`, `{d}/non_consortium_sequences.csv` | — |
| merge_excludes | (inline in `run.py`) | `exclude.txt`, `{d}/exclude_auto.txt` | `{d}/exclude_merged.txt` | — |
| mafft *(optional)* | `mafft` | `aligned.fasta` | `{d}/{p}.mafft.fasta` | `mafft` |
| clean | `clean_up_pathoplexus.py` | `aligned.fasta`, `metadata.tsv`, `aliases.csv`, `mask_sites.txt` (optional) | `{d}/{p}.cleaned.aln` | — |
| qc_basic | `genome_qc_report.py` | `{d}/{p}.cleaned.aln`, `metadata.tsv` | `{d}/qc_report_basic.txt` | — |
| pre_filter | (inline in `run_pipeline.py`) | `{d}/{p}.cleaned.aln`, `{d}/exclude_merged.txt`, exclude_quality (optional) | `{d}/{p}.aln` | — |
| iqtree1 | `iqtree3` | `{d}/{p}.aln` | `{d}/{p}.aln.treefile`, `{d}/{p}.aln.state`, `{d}/{p}.aln.iqtree`, `{d}/best_model.txt` | `iqtree3` |
| raccoon | `raccoon` | `{d}/{p}.aln.treefile`, `{d}/{p}.aln`, `{d}/{p}.aln.state` | `{d}/{p}.aln.phylo_flags.csv`, `{d}/{p}.aln.branch_snps.reconstruction.csv` | `raccoon` |
| adar_mask | `mask_adar_sites.py` | `{d}/{p}.aln`, `phylo_flags.csv`, `.treefile`, `.state` | `{d}/{p}.ADAR_masked.aln`, `{d}/{p}.ADAR_masked.adar_report.txt` | — |
| qc_phylo | `phylo_qc_report.py` | `phylo_flags.csv`, `branch_snps.reconstruction.csv`, `metadata.tsv`, `qc_report_basic.txt` | `{d}/qc_report_phylo.txt` | — |
| iqtree2 | `iqtree3` | `{d}/{p}.ADAR_masked.aln`, `best_model.txt` | `{d}/{p}.ADAR_masked.aln.treefile` | `iqtree3` |
| temporal_root | `temporal_root.py` | `{d}/{p}.ADAR_masked.aln.treefile` | `{d}/{p}.ADAR_masked.bestFitRoot.tree`, `{d}/{p}.ADAR_masked.outliers.txt`, `{d}/{p}.ADAR_masked.rtt_plot.pdf`, `{d}/{p}.ADAR_masked.rtt_summary.txt` | — |
| filter_outliers | (inline in `run_pipeline.py`) | `{d}/{p}.ADAR_masked.aln`, `{d}/outliers.txt` | `{d}/{p}.ADAR_masked.final.aln` | — |
| iqtree3 | `iqtree3` | `{d}/{p}.ADAR_masked.final.aln`, `best_model.txt` | `{d}/{p}.ADAR_masked.final.aln.treefile` | `iqtree3` |
| final_root | `temporal_root.py` | `{d}/{p}.ADAR_masked.final.aln.treefile` | `{d}/{p}.ADAR_masked.final.bestFitRoot.tree`, `{d}/{p}.ADAR_masked.final.rtt_plot.pdf`, `{d}/{p}.ADAR_masked.final.rtt_summary.txt` | — |
| health_checks | `health_checks.py` | `{d}/` (current run dir), previous run dir (optional) | `{d}/health_check_report.txt` | — |
| compare_runs | `compare_runs.py` | `{d}/` (current), previous run dir | `{d}/run_comparison.txt` | — |
| beast_xmls | `generate_beast_xml.py` | `{d}/{p}.ADAR_masked.final.aln`, template XMLs | `{d}/*.xml` | — |

**Step dependency chain:**

```
download → metadata_qc → merge_excludes
                                ↓
         mafft (optional) → clean → qc_basic
                                ↓
                          pre_filter → iqtree1 → raccoon → adar_mask → qc_phylo
                                                              ↓
                                                          iqtree2 → temporal_root → filter_outliers → iqtree3 → final_root
                                                                                                                     ↓
                                                                                                              beast_xmls
                                                                                                              health_checks
                                                                                                              compare_runs
```

**Notes for Nextflow:**
- All Python scripts exit with code 1 on failure.
- `iqtree3` must be on `PATH`. Three IQ-TREE calls: `iqtree1` uses `--ancestral` + either `MFP` or a fixed model; `iqtree2` uses the saved model; `iqtree3` adds `-abayes`.
- `raccoon` must be on `PATH` (installed via pip in the conda env).
- `qc_basic` and `qc_phylo` are reporting-only steps with no downstream dependencies — they can run in parallel with the main chain.
- `health_checks`, `compare_runs`, and `beast_xmls` are independent of each other and can run in parallel after `final_root`.
- `merge_excludes` is an inline function in `run.py` — it deduplicates and concatenates text files. Will need to be extracted or reimplemented for Nextflow.
- `pre_filter` and `filter_outliers` are inline functions in `run_pipeline.py` — they filter a FASTA by accession IDs in exclude files. Will need to be extracted or reimplemented for Nextflow.

## Environment notes

- **raccoon** requires `baltic<1` (pinned in `environment.yml`)
- **IQ-TREE** must be installed separately (not available via conda for all platforms)
