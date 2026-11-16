// iqtree1: ML tree + ModelFinder (or fixed --model) + ancestral reconstruction
process iqtree_ancestral {

  input:
    path alignment

  output:
    tuple path(alignment), path("${alignment}.treefile"), path("${alignment}.state"), emit: tree
    path "best_model.txt", emit: model
    path "${alignment}.iqtree"

  script:
    model = params.model ?: "MFP"
    """
    # threads: the requested number, limited to what this machine or container allows
    nt=\$(effective_cpus.sh ${task.cpus})
    if [ "\$nt" -gt 1 ]; then threads="-nt AUTO -ntmax \$nt"; else threads="-nt 1"; fi
    echo "IQ-TREE threads: \$threads (requested up to ${task.cpus})"
    iqtree3 \$threads -s ${alignment} -m ${model} -czb --ancestral -redo -seed ${params.seed}
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

  input:
    path alignment
    val model

  output:
    path "${alignment}.treefile", emit: tree
    path "${alignment}.iqtree"

  script:
    args = task.ext.args ?: ""
    """
    # threads: the requested number, limited to what this machine or container allows
    nt=\$(effective_cpus.sh ${task.cpus})
    if [ "\$nt" -gt 1 ]; then threads="-nt AUTO -ntmax \$nt"; else threads="-nt 1"; fi
    echo "IQ-TREE threads: \$threads (requested up to ${task.cpus})"
    iqtree3 \$threads -s ${alignment} -m ${model} -czb -redo -seed ${params.seed} ${args}
    """
}
