process beast_xml {

  input:
    path alignment
    path templates
    val prefix

  output:
    path "*.xml"

  script:
    extra = ""
    if ( params.template_egc )
        extra += " --template-egc ${file(params.template_egc).name}"
    if ( params.template_sg )
        extra += " --template-sg ${file(params.template_sg).name} --anchor-date ${params.anchor_date}"
    if ( params.template_lgc )
        extra += " --template-lgc ${file(params.template_lgc).name}"
    if ( params.chain_length )
        extra += " --chain-length ${params.chain_length}"
    if ( params.log_every )
        extra += " --log-every ${params.log_every}"
    """
    n=\$(grep -c '>' ${alignment})
    generate_beast_xml.py -a ${alignment} -o . -p ${prefix}.n\$n ${extra}
    """
}
