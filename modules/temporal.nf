// Imported as temporal_root (with --iterate, from ext.args) and final_root
process temporal_root {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path tree
    val prefix

  output:
    path "${prefix}.outliers.txt", emit: outliers
    path "${prefix}.*", emit: all

  script:
    args = task.ext.args ?: ""
    """
    export MPLCONFIGDIR=\$PWD/.mpl
    temporal_root.py -t ${tree} --plot -o . --prefix ${prefix} ${args}
    """
}

process select_outliers {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path detected_outliers

  output:
    path "outliers.txt"

  script:
    """
    select_outliers.py -i ${detected_outliers} -t ${params.outlier_removal_threshold} -o outliers.txt
    """
}
