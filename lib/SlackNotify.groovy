import groovy.json.JsonOutput
import nextflow.Nextflow

/**
 * Posts a short summary to a Slack incoming webhook (--slack_webhook, or the SLACK_WEBHOOK_URL
 * environment variable) when the run finishes or fails. A failure to notify never fails the run.
 * The webhook URL is never printed or logged.
 */
class SlackNotify {

    static void send(workflow, Map params, log) {
        def url = params.slack_webhook
        if ( !url || url == false )
            return
        try {
            def lines = []
            lines << "BDBV pipeline ${workflow.success ? 'finished' : 'FAILED'}: ${params.outdir} (${params.prefix})".toString()
            def health = Nextflow.file("${params.outdir}/11_health_checks/health_check_report.txt")
            if ( workflow.success && health.exists() ) {
                def summary = health.text.readLines().findAll { it.trim() && !it.startsWith('=') }.take(20)
                lines << '```' + summary.join('\n') + '```'
            } else if ( !workflow.success ) {
                lines << "See ${params.outdir}/${params.prefix}.run_log.txt for the output of each step".toString()
            }
            def conn = new URL(url.toString()).openConnection()
            conn.setRequestMethod('POST')
            conn.setConnectTimeout(10000)
            conn.setReadTimeout(10000)
            conn.setDoOutput(true)
            conn.setRequestProperty('Content-Type', 'application/json')
            conn.outputStream.withWriter('UTF-8') { it << JsonOutput.toJson([text: lines.join('\n')]) }
            def code = conn.responseCode
            if ( code >= 200 && code < 300 )
                log.info "Slack notification sent"
            else
                log.warn "Slack notification failed (HTTP ${code})"
        } catch ( Exception e ) {
            // do not include the URL: java exceptions can echo it
            log.warn "Slack notification failed (${e.class.simpleName})"
        }
    }
}
