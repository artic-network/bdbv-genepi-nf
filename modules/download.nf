process download {

  output:
    path "${params.organism}_aligned-nuc_${params.run_name}.fasta", emit: fasta
    path "${params.organism}_metadata_${params.run_name}.tsv", emit: metadata

  script:
    extra = ""
    if ( params.max_date )
        extra += " --max-date ${params.max_date}"
    // free text: quote it for the shell (a ' inside becomes '\'')
    if ( params.pathoplexus_filters )
        extra += " --filters '" + params.pathoplexus_filters.toString().replace("'", "'\\''") + "'"
    """
    download_pathoplexus.py -o download --organism ${params.organism} --scope ${params.scope}${extra}
    mv download/aligned.fasta ${params.organism}_aligned-nuc_${params.run_name}.fasta
    mv download/metadata.tsv ${params.organism}_metadata_${params.run_name}.tsv
    """
}
