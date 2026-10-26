// iqtree1: ML tree + ModelFinder (or fixed --model) + ancestral reconstruction
process iqtree_ancestral {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path alignment

  output:
    tuple path(alignment), path("${alignment}.treefile"), path("${alignment}.state"), emit: tree
    path "best_model.txt", emit: model
    path "${alignment}.iqtree"

  script:
    model = params.model ?: "MFP"
    """
    iqtree3 -nt ${task.cpus} -s ${alignment} -m ${model} -czb --ancestral -redo -seed ${params.seed}
    if [ "${model}" = "MFP" ]; then
        sed -n 's/^Best-fit model according to BIC:[[:space:]]*//p' ${alignment}.iqtree > best_model.txt
        [ -s best_model.txt ] || { echo "ERROR: could not find best-fit model in ${alignment}.iqtree" >&2; exit 1; }
    else
        echo "${model}" > best_model.txt
    fi
    """
}

// Imported as iqtree_adar (iqtree2) and iqtree_final (iqtree3); extra flags come from ext.args
process iqtree {

  publishDir "${params.outdir}", mode: params.publish_dir_mode

  input:
    path alignment
    val model

  output:
    path "${alignment}.treefile", emit: tree
    path "${alignment}.iqtree"

  script:
    args = task.ext.args ?: ""
    """
    iqtree3 -nt ${task.cpus} -s ${alignment} -m ${model} -czb -redo -seed ${params.seed} ${args}
    """
}
