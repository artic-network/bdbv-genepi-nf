// BDBV phylogenetic pipeline: Pathoplexus download -> cleaned, ADAR-handled, outlier-filtered
// ML tree -> health checks. Nextflow port of bin/run.py and bin/run_pipeline.py.

include { download } from './modules/download'
include { mafft; clean; qc_basic; metadata_quality_report; merge_exclude } from './modules/prepare'
include { filter_exclude as pre_filter; filter_exclude as filter_outliers } from './modules/prepare'
include { iqtree_ancestral } from './modules/iqtree'
include { iqtree as iqtree_adar; iqtree as iqtree_final } from './modules/iqtree'
include { raccoon; adar_mask; qc_phylo } from './modules/adar'
include { temporal_root; temporal_root as final_root; select_outliers } from './modules/temporal'
include { health_checks; compare_runs } from './modules/health'

workflow {
  adar_tag = params.adar_mode == 'mask' ? 'ADAR_masked' : 'ADAR_excluded'
  adar_prefix = "${params.prefix}.${adar_tag}"

  // Inputs: download from Pathoplexus unless a FASTA and metadata are given
  if ( params.fasta ) {
      if ( !params.metadata )
          error "--metadata is required when --fasta is given"
      fasta_ch = channel.value(file("${params.fasta}", checkIfExists:true))
      metadata_ch = channel.value(file("${params.metadata}", checkIfExists:true))
  } else {
      download()
      fasta_ch = download.out.fasta
      metadata_ch = download.out.metadata
  }

  aliases_ch = channel.value(file("${params.aliases_file}", checkIfExists:true))
  manual_exclude_ch = channel.fromPath([params.exclude_file].findAll { it }, checkIfExists:true)

  // INRB linelist: duplicate resolution, mismatch checks and automatic exclusions
  auto_exclude_ch = channel.empty()
  if ( params.linelist ) {
      linelist_ch = channel.value(file("${params.linelist}", checkIfExists:true))
      hz_ch = channel.value(params.hz_lookup ? file("${params.hz_lookup}", checkIfExists:true) : [])
      metadata_quality_report(metadata_ch, fasta_ch, aliases_ch, linelist_ch, hz_ch)
      auto_exclude_ch = metadata_quality_report.out.exclude_auto
  }
  pending_exclude_ch = channel.fromPath([params.exclude_pending].findAll { it }, checkIfExists:true)

  // pre_filter removes manual, automatic and pending exclusions
  exclude_ch = manual_exclude_ch.mix(auto_exclude_ch, pending_exclude_ch)
      .collect()
      .ifEmpty([])
  // The merged list (manual + automatic) is kept as a record and read by compare_runs
  merged_exclude_ch = merge_exclude(manual_exclude_ch.mix(auto_exclude_ch).collect())

  // Clean, QC and exclusions
  if ( params.mafft ) {
      fasta_ch = mafft(fasta_ch)
  }
  mask_sites_ch = channel.value(params.mask_sites_file ? file("${params.mask_sites_file}", checkIfExists:true) : [])
  clean(fasta_ch, metadata_ch, aliases_ch, mask_sites_ch)
  qc_basic(clean.out.aln, metadata_ch)
  pre_filter(clean.out.aln, exclude_ch, params.prefix)

  // Initial tree, ADAR identification and handling
  iqtree_ancestral(pre_filter.out)
  model_ch = iqtree_ancestral.out.model.map { it.text.trim() }
  raccoon(iqtree_ancestral.out.tree)
  adar_mask(iqtree_ancestral.out.tree, raccoon.out.flags, adar_prefix)
  qc_phylo(raccoon.out.flags, metadata_ch, qc_basic.out)

  iqtree_adar(adar_mask.out.aln, model_ch)
  temporal_root(iqtree_adar.out.tree, adar_prefix)

  // Review point 1, outliers: auto-remove at a z threshold, use a reviewed outlier file, or stop here
  outliers_ch = null
  if ( params.outlier_removal_threshold != null ) {
      outliers_ch = select_outliers(temporal_root.out.outliers)
  } else if ( params.outlier_file ) {
      outliers_ch = channel.value(file("${params.outlier_file}", checkIfExists:true))
  }

  if ( outliers_ch != null ) {
      filter_outliers(adar_mask.out.aln, outliers_ch, "${adar_prefix}.final")
      iqtree_final(filter_outliers.out, model_ch)
      final_root(iqtree_final.out.tree, "${adar_prefix}.final")

      // Health checks against a previous run (update runs only)
      if ( params.previous_data ) {
          current_files = temporal_root.out.all
              .mix(final_root.out.all, adar_mask.out.aln, adar_mask.out.report, raccoon.out.flags.flatten())
              .collect()
          previous_ch = channel.value(file("${params.previous_data}", type: "dir", checkIfExists:true))
          pending_ch = channel.value(params.exclude_pending ? file("${params.exclude_pending}", checkIfExists:true) : [])
          health_checks(current_files, metadata_ch, previous_ch, pending_ch)

          // Run comparison against the previous run
          comparison_ch = clean.out.aln
              .mix(pre_filter.out, adar_mask.out.aln, adar_mask.out.report, filter_outliers.out, final_root.out.all,
                   temporal_root.out.all, merged_exclude_ch)
              .mix(outliers_ch.flatten())
          if ( params.linelist )
              comparison_ch = comparison_ch.mix(metadata_quality_report.out.csv)
          compare_runs(comparison_ch.collect(), metadata_ch, previous_ch)
      }
  }
}

// Message for the outlier review point. It is printed on completion so it is the last thing on screen.
workflow.onComplete {
  RunLog.write(workflow, params)
  SlackNotify.send(workflow, params, log)
  if ( !workflow.success )
      return
  def adar_tag = params.adar_mode == 'mask' ? 'ADAR_masked' : 'ADAR_excluded'
  def outliers = "${params.outdir}/08_temporal_root/${params.prefix}.${adar_tag}.outliers.txt"
  if ( params.outlier_removal_threshold == null && !params.outlier_file ) {
      log.info """
      ==========================================================================
      Stopped after temporal rooting for outlier review
      (no --outlier_removal_threshold or --outlier_file given).

      1. Review the RTT plot and detected outliers in ${params.outdir}/08_temporal_root/
      2. Edit a copy of: ${outliers}
      3. Resume, passing your edited file:

         ${workflow.commandLine} -resume --outlier_file <edited_outliers_file>
      ==========================================================================
      """.stripIndent()
  }
}

// The run log is also written when the pipeline fails, so the output of the failed step is kept.
workflow.onError {
  RunLog.write(workflow, params)
}
