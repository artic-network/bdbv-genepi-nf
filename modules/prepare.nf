process mafft {

  input:
    path fasta

  output:
    path "${params.prefix}.mafft.fasta"

  script:
    """
    # An input that is already aligned (for example from Pathoplexus) is not re-aligned: MAFFT would
    # change the number of columns, and the trim and mask positions refer to the original alignment
    # (exit 0 = aligned, 1 = not aligned; anything else means the check itself failed, so stop)
    status=0
    check_alignment.py --quiet ${fasta} || status=\$?
    if [ "\$status" -eq 0 ]; then
        echo "Input is already aligned (all sequences the same length): skipping MAFFT and using it as it is"
        cp -L ${fasta} ${params.prefix}.mafft.fasta
    elif [ "\$status" -eq 1 ]; then
        mafft --auto --thread -1 ${fasta} > ${params.prefix}.mafft.fasta
    else
        echo "ERROR: could not check whether ${fasta} is aligned (check_alignment.py exited with status \$status)" >&2
        exit 1
    fi
    """
}

process clean {

  input:
    path fasta
    path metadata
    path aliases
    path mask_sites_file

  output:
    path "${params.prefix}.cleaned.aln", emit: aln
    path "${params.prefix}.cleaned.log"
    path "missing_aliases.txt", optional: true
    path "${params.prefix}.cleaned.flagged_early_dates.txt", optional: true

  script:
    extra = ""
    if ( params.trim_start != null )
        extra += " --trim_start ${params.trim_start}"
    if ( params.trim_end != null )
        extra += " --trim_end ${params.trim_end}"
    if ( params.min_date )
        extra += " --min-date ${params.min_date}"
    if ( params.mask_sites )
        extra += " --mask_sites ${params.mask_sites}"
    if ( params.mask_sites_file )
        extra += " --mask_sites_file ${mask_sites_file}"
    hint = params.mafft ? "" : "--mafft-hint"
    """
    # fail early, rather than when IQ-TREE reads the alignment, if the input is not aligned
    check_alignment.py ${fasta} ${hint}
    clean_up_pathoplexus.py -f ${fasta} -m ${metadata} -a ${aliases} -p ${params.prefix}.cleaned --missing-aliases-out missing_aliases.txt ${extra} \
        | tee ${params.prefix}.cleaned.log
    """
}

process qc_basic {

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

// Checks the Pathoplexus metadata (duplicates, locations, dates, sequence quality) and lists accessions
// to exclude automatically. With the optional INRB linelist (--linelist) duplicates are also resolved
// against it and date/location mismatches against it are reported.
process metadata_quality_report {

  input:
    path metadata
    path fasta
    path aliases
    path linelist
    path hz_lookup

  output:
    path "metadata_quality_report.txt", emit: report
    path "metadata_quality_report.csv", emit: csv
    path "exclude_auto.txt", emit: exclude_auto
    path "non_consortium_sequences.csv", optional: true

  script:
    hz = params.hz_lookup ? "-s ${hz_lookup}" : ""
    ll = params.linelist ? "-l ${linelist}" : ""
    """
    metadata_quality_report.py -m ${metadata} -f ${fasta} -a ${aliases} ${ll} ${hz} \
        -o metadata_quality_report.csv --exclude-out exclude_auto.txt \
        --report-txt metadata_quality_report.txt --non-consortium-out non_consortium_sequences.csv
    # the script only writes the exclude list when there is something to exclude
    [ -f exclude_auto.txt ] || echo "# nothing to exclude automatically" > exclude_auto.txt
    """
}

// Manual plus automatic exclusions as one deduplicated list (compare_runs reads it from each run)
process merge_exclude {

  input:
    path exclude_files

  output:
    path "exclude_merged.txt"

  script:
    """
    echo "# Merged exclude list (manual + auto)" > exclude_merged.txt
    cat ${exclude_files} | grep -v '^[[:space:]]*#' | awk 'NF {print \$1}' | awk '!seen[\$0]++' >> exclude_merged.txt
    """
}
