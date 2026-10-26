process mafft {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path fasta

  output:
    path "${params.prefix}.mafft.fasta"

  script:
    """
    mafft --auto --thread ${task.cpus} ${fasta} > ${params.prefix}.mafft.fasta
    """
}

process clean {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path fasta
    path metadata
    path aliases

  output:
    path "${params.prefix}.cleaned.aln", emit: aln
    path "${params.prefix}.cleaned.log"
    path "${params.prefix}.cleaned.flagged_early_dates.txt", optional: true

  script:
    extra = ""
    if ( params.trim_start != null )
        extra += " --trim_start ${params.trim_start}"
    if ( params.trim_end != null )
        extra += " --trim_end ${params.trim_end}"
    if ( params.min_date )
        extra += " --min-date ${params.min_date}"
    """
    clean_up_pathoplexus.py -f ${fasta} -m ${metadata} -a ${aliases} -p ${params.prefix}.cleaned ${extra} \
        | tee ${params.prefix}.cleaned.log
    """
}

process qc_basic {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path alignment
    path metadata

  output:
    path "qc_report_basic.txt"

  script:
    """
    genome_qc_report.py -a ${alignment} -m ${metadata} -o qc_report_basic.txt
    """
}

// Used as pre_filter (exclusion lists) and filter_outliers (temporal outliers)
process filter_exclude {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path alignment
    path exclude_files
    val prefix

  output:
    path "${prefix}.aln"

  script:
    """
    filter_exclude.py --input ${alignment} --output ${prefix}.aln --exclude-files ${exclude_files}
    """
}
