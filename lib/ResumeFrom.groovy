import nextflow.Nextflow
import nextflow.exception.AbortOperationException

/**
 * Finds the outputs of an earlier run that stopped for outlier review (--resume_from), so a run can
 * continue from filter_outliers without repeating the download, the trees or the ADAR steps.
 * The files are looked up in the numbered output folders by their names, which contain the run's
 * --prefix and --adar_mode, so those must match the earlier run.
 */
class ResumeFrom {

    static def one(String dir, String relative) {
        def f = Nextflow.file("${dir}/${relative}")
        if ( !f.exists() ) {
            def folder = relative.contains('/') ? relative.substring(0, relative.lastIndexOf('/')) : ''
            def found = Nextflow.files("${dir}/${folder}/*").collect { it.name }.sort()
            throw new AbortOperationException(
                "--resume_from: cannot find ${dir}/${relative}\n" +
                "  --resume_from must be the output folder of an earlier run that reached temporal rooting, " +
                "and --prefix and --adar_mode must match that run.\n" +
                "  Files in ${dir}/${folder}: ${found ? found.join(', ') : '(none, or the folder is missing)'}")
        }
        return f
    }

    static List many(String dir, String pattern) {
        return Nextflow.files("${dir}/${pattern}").sort()
    }

    /** First file matching a pattern, or null (used for the optional metadata copy). */
    static def first(String dir, String pattern) {
        def l = many(dir, pattern)
        return l ? l[0] : null
    }
}
