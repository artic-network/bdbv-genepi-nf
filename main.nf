process mafft {

  publishDir "${params.outdir}", mode: 'copy'

  input:
    path fasta

  output:
    path "${fasta.baseName}.aln"

  script:
    """
    mafft --auto --thread -1 ${fasta} > ${fasta.baseName}.aln
    """
}

process clean {
  
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path fasta
    path metadata
    
  output:
    path "${fasta.baseName}.cleaned.aln"
    path "${fasta.baseName}.cleaned.log"
    
  script:
    extra = ""
    if ( params.min_date )
        extra += " --min-date ${params.min_date}"
    """
    python clean.py \
        --fasta ${fasta} \
        --metadata ${metadata} \
        -a ${params.aliases_file} \
        --prefix ${fasta.baseName}.cleaned \
        ${extra} > ${fasta.baseName}.cleaned.log
    """
}

process qc_basic {

  publishDir "${params.outdir}", mode: 'copy'

  input:
    path alignment
    path metadata

  output:
    path "qc_basic_report.txt"

  script:
    """
    python qc_basic.py -a ${alignment} -m ${metadata} --output qc_basic_report.txt
    """
}

process filter_exclude {

  publishDir "${params.outdir}", mode: 'copy'

  input:
    path alignment
    path exclude_identity
    path exclude_quality

  output:
    path "${alignment.baseName}.filtered.aln"

  script:
    """
    python filter_exclude.py \
        --input ${alignment} \
        --output ${alignment.baseName}.filtered.aln \
        --exclude-files ${exclude_identity} ${exclude_quality}
    """
}


workflow {
    fasta_ch = channel.of(file("${params.fasta}", type: "file", checkIfExists:true))
    metadata_ch = channel.of(file("${params.metadata}", type: "file", checkIfExists:true))
    
    if ( params.mafft ) {
        mafft_ch = mafft(fasta_ch)
        clean(mafft_ch, metadata_ch)
    } else {
        cleaned_ch = clean(fasta_ch, metadata_ch)
    }

    qc_basic(cleaned_ch, metadata_ch)

    exclude_identity = file("${params.exclude_file}", type: "file", checkIfExists: true)
    exclude_quality = file("${params.exclude_quality}", type: "file", checkIfExists: true)
    filter_exclude(cleaned_ch, exclude_identity, exclude_quality) 
}
