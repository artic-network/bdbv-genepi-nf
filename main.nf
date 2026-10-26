// BDBV phylogenetic pipeline: Pathoplexus download -> cleaned, ADAR-handled, outlier-filtered
// ML tree -> health checks -> BEAST XMLs. Nextflow port of bin/run_pipeline.py and bin/run_update.py.

include { download } from './modules/download'
include { mafft; clean; qc_basic } from './modules/prepare'
include { filter_exclude as pre_filter; filter_exclude as filter_outliers } from './modules/prepare'
include { iqtree_ancestral } from './modules/iqtree'
include { iqtree as iqtree_adar; iqtree as iqtree_final } from './modules/iqtree'
include { raccoon; adar_mask; qc_phylo } from './modules/adar'
include { temporal_root; temporal_root as final_root; select_outliers } from './modules/temporal'
include { health_checks } from './modules/health'
include { beast_xml } from './modules/beast'

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
  exclude_ch = channel.fromPath([params.exclude_file, params.exclude_quality, params.exclude_pending].findAll { it }, checkIfExists:true)
      .collect()
      .ifEmpty([])

  // Clean, QC and exclusions
  if ( params.mafft ) {
      fasta_ch = mafft(fasta_ch)
  }
  clean(fasta_ch, metadata_ch, aliases_ch)
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
  } else {
      log.info "Stopping after temporal rooting for outlier review (no --outlier_removal_threshold or --outlier_file)."
      log.info "Edit a copy of ${params.outdir}/${adar_prefix}.outliers.txt, then rerun with -resume --outlier_file <file>"
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
      }

      // Review point 2, health checks: optionally stop before BEAST XMLs
      if ( params.previous_data && params.stop_after_health_checks ) {
          log.info "Stopping after health checks for review. Rerun with -resume and without --stop_after_health_checks to generate BEAST XMLs."
      } else if ( params.beast_xml ) {
          templates_ch = channel.fromPath([params.template_egc, params.template_sg, params.template_lgc].findAll { it }, checkIfExists:true)
              .collect()
          beast_xml(filter_outliers.out, templates_ch, adar_prefix)
      }
  }
}
