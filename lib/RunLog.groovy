import nextflow.Nextflow

/**
 * Builds one run log for the whole pipeline: a header (command, parameters, outcome) followed by the
 * screen output (stdout and stderr) of every task, in order. It is written whether the run succeeds
 * or fails, so the output of the step that failed is included. Tasks are found from the trace file
 * (see the trace block in nextflow.config).
 */
class RunLog {

    static void write(workflow, Map params) {
        def outdir = params.outdir.toString()
        def trace = Nextflow.file("${outdir}/execution/trace.txt")
        def out = new StringBuilder()

        def tasks = []
        if ( trace.exists() ) {
            def lines = trace.text.readLines().findAll { it.trim() }
            def cols = lines ? lines[0].split('\t', -1).toList() : []
            tasks = lines.drop(1).collect { l -> [cols, l.split('\t', -1).toList()].transpose().collectEntries { [(it[0]): it[1]] } }
            tasks = tasks.sort { it.task_id as Integer }
        }
        def failed = tasks.findAll { it.status == 'FAILED' }

        out << "PIPELINE RUN LOG\n"
        out << "Started:   ${workflow.start}\n"
        out << "Finished:  ${workflow.complete ?: new Date()}\n"
        out << "Status:    ${workflow.success ? 'SUCCESS' : 'FAILED'}\n"
        if ( failed )
            out << "Failed:    ${failed.collect { "${it.name} (exit ${it.exit})" }.join(', ')}\n"
        out << "Command:   ${workflow.commandLine}\n"
        out << "Output:    ${outdir}\n"
        out << "Prefix:    ${params.prefix}\n"
        out << "Pipeline:  ${workflow.manifest?.name ?: ''} ${workflow.revision ?: ''} (Nextflow ${workflow.nextflow.version})\n"
        out << "Container: ${workflow.container ?: '(none)'}\n"
        out << "\nParameters:\n"
        params.sort().each { k, v ->
            if ( v != null && !(v instanceof Map) && !(k in ['help', 'version', 'email', 'email_on_fail', 'plaintext_email',
                    'monochrome_logs', 'hook_url', 'validate_params', 'show_hidden_params', 'schema_ignore_params']) )
                out << "  ${k}: ${v}\n"
        }
        out << "\n"

        if ( trace.exists() ) {
            tasks.each { t ->
                out << "── ${t.name} ── [${t.submit}] ${t.status}, exit ${t.exit}\n"
                def logf = t.workdir ? Nextflow.file("${t.workdir}/.command.log") : null
                if ( logf && logf.exists() && logf.size() > 0 )
                    out << logf.text.trim() << "\n"
                else
                    out << "  (no output)\n"
                out << "\n"
            }
        } else {
            out << "(no trace file found, so no task output could be collected)\n"
        }

        // Both onComplete and onError run on failure: replace the file rather than overwrite in place
        def logFile = Nextflow.file("${outdir}/${params.prefix}.run_log.txt")
        logFile.delete()
        logFile.text = out.toString()
    }
}
