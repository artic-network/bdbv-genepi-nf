process download {

  output:
    path "${params.organism}_aligned-nuc_${params.run_name}.fasta", emit: fasta
    path "${params.organism}_metadata_${params.run_name}.tsv", emit: metadata

  script:
    extra = ""
    if ( params.max_date )
        extra += " --max-date ${params.max_date}"
    """
    download_pathoplexus.py -o download --organism ${params.organism} --scope ${params.scope}${extra}
    mv download/aligned.fasta ${params.organism}_aligned-nuc_${params.run_name}.fasta
    mv download/metadata.tsv ${params.organism}_metadata_${params.run_name}.tsv
    """
}
