process raccoon {

  input:
    tuple path(alignment), path(treefile), path(statefile)

  output:
    tuple path("raccoon_output/*.phylo_flags.csv"), path("raccoon_output/*.branch_snps.reconstruction.csv"), emit: flags
    path "raccoon_output/**"

  script:
    """
    export XDG_CACHE_HOME=\$PWD/.cache
    raccoon tree-qc -t ${treefile} --alignment ${alignment} --asr-state ${statefile} --run-adar -d raccoon_output
    """
}

process adar_mask {

  input:
    tuple path(alignment), path(treefile), path(statefile)
    tuple path(phylo_flags), path(branch_snps)
    val prefix

  output:
    path "${prefix}.aln", emit: aln
    path "${prefix}.adar_report.txt", emit: report
    path "${prefix}.adar_exclusions.txt", optional: true

  script:
    """
    mask_adar_sites.py -a ${alignment} -r ${phylo_flags} -t ${treefile} -s ${statefile} --mode ${params.adar_mode} -p ${prefix}
    """
}

process qc_phylo {

  input:
    tuple path(phylo_flags), path(branch_snps)
    path metadata
    path qc_basic_report

  output:
    path "qc_report_phylo.txt"

  script:
    """
    phylo_qc_report.py -f ${phylo_flags} -b ${branch_snps} -m ${metadata} --basic-qc ${qc_basic_report} -o qc_report_phylo.txt
    """
}
