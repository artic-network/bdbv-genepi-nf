// Compares this run with a previous run directory. Writes suggested additions to the
// pending exclusion list and new retractions to outdir instead of editing files in data/.
process health_checks {

  input:
    path current_files, stageAs: "current/*"
    path metadata, stageAs: "current/bdbv_metadata_current.tsv"
    path previous, stageAs: "previous"
    path exclude_pending

  output:
    path "health_check_report.txt", emit: report
    path "exclude_pending.additions.txt"
    path "retractions.additions.txt", optional: true

  script:
    full = params.full_health_check ? "--full" : ""
    """
    health_checks.py --current current --previous previous --prefix ${params.prefix} ${full} -o health_check_report.txt
    health_check_pending.py --current current --previous previous --prefix ${params.prefix} ${full} \
        --exclude-pending ${exclude_pending} -o exclude_pending.additions.txt
    if [ -f retractions.txt ]; then
        sed 's/\tcurrent\t/\t${params.run_name}\t/' retractions.txt > retractions.additions.txt
    fi
    """
}
