import nextflow.Nextflow

/**
 * Builds one run log for the whole pipeline: a header (command, parameters, outcome) followed by the
 * screen output (stdout and stderr) of every task, in order. It is written whether the run succeeds
 * or fails, so the output of the step that failed is included. Tasks are found from the trace file
 * (see the trace block in nextflow.config).
 */
class RunLog {

    /** Never write the Slack webhook URL (a secret) to the log. */
    static String redact(String text) {
        return text?.replaceAll(/(--slack_webhook[ =])\S+/, '$1<redacted>')
    }

    /** Exit status of a task in words: out of memory and out of time are the common causes of a stop. */
    static String describe(Map t) {
        def exit = t.exit
        if ( exit == '137' )
            return 'exit 137: killed, usually because it ran out of memory'
        if ( exit == '143' )
            return 'exit 143: terminated'
        if ( exit in [null, '', '-'] && t.status == 'FAILED' )
            return 'no exit status: stopped before finishing, usually because it exceeded its time limit' + (t.time ? " (${t.time})" : '')
        return "exit ${exit}"
    }

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
            out << "${workflow.success ? 'Retried:   ' : 'Failed:    '}${failed.collect { "${it.name} (${describe(it)})" }.join(', ')}\n"
        out << "Command:   ${redact(workflow.commandLine)}\n"
        out << "Output:    ${outdir}\n"
        out << "Prefix:    ${params.prefix}\n"
        out << "Pipeline:  ${workflow.manifest?.name ?: ''} ${workflow.revision ?: ''} (Nextflow ${workflow.nextflow.version})\n"
        out << "Container: ${workflow.container ?: '(none)'}\n"
        out << "Host:      ${System.getProperty('os.name')} ${System.getProperty('os.version')} ${System.getProperty('os.arch')}, " +
               "${Runtime.runtime.availableProcessors()} CPUs, Java ${System.getProperty('java.version')}\n"
        out << "\nParameters:\n"
        params.sort().each { k, v ->
            if ( v != null && !(v instanceof Map) && !(k in ['help', 'version', 'email', 'email_on_fail', 'plaintext_email',
                    'monochrome_logs', 'hook_url', 'validate_params', 'show_hidden_params', 'schema_ignore_params', 'slack_webhook']) )
                out << "  ${k}: ${v}\n"
        }
        out << "\n"

        // On failure, include Nextflow's own error report: the message shown on screen, with the failed
        // command and its error output. It does not depend on the trace file, which may not have been
        // written yet (or may not list the failed task) when the run log is built.
        if ( !workflow.success && workflow.errorReport ) {
            out << "── Nextflow error report ──\n"
            out << redact(workflow.errorReport.toString().trim()) << "\n\n"
        }

        if ( trace.exists() ) {
            tasks.each { t ->
                out << "── ${t.name} ── [${t.submit}] ${t.status}, ${describe(t)}\n"
                out << "  attempt ${t.attempt ?: '?'}; requested ${t.memory ?: '?'} memory, ${t.time ?: '?'} time; " +
                       "ran ${t.realtime ?: '?'}, peak memory ${t.peak_rss ?: '?'}\n"
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

        // Nextflow's own log (Docker errors, executor messages) is the last place a cause can hide
        if ( !workflow.success ) {
            def logPath = (workflow.commandLine =~ /-log\s+('([^']+)'|(\S+))/).with { it.find() ? (it.group(2) ?: it.group(3)) : null }
            def nfLog = Nextflow.file(logPath ?: "${workflow.launchDir}/.nextflow.log")
            if ( nfLog.exists() ) {
                out << "── Last lines of ${nfLog.name} ──\n"
                out << redact(nfLog.text.readLines().takeRight(60).join("\n")) << "\n\n"
            }
        }

        // Both onComplete and onError run on failure: replace the file rather than overwrite in place
        def logFile = Nextflow.file("${outdir}/${params.prefix}.run_log.txt")
        logFile.delete()
        logFile.text = out.toString()
    }
}
