/*
 * modules/checkv.nf
 * CHECKV — genome completeness and quality assessment.
 * Reused from wvdb_build; same process definition.
 */

process CHECKV {
    label 'cpu_medium'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path input_fasta
    path checkvdb

    output:
    path "checkv_out/quality_summary.tsv",   emit: quality_summary
    path "checkv_out/completeness.tsv",       emit: completeness
    path "checkv_out/contamination.tsv",      emit: contamination
    path "checkv_out/",                        emit: checkv_dir

    script:
    """
    checkv end_to_end \\
        -d "${checkvdb}" \\
        -t ${task.cpus} \\
        --remove_tmp \\
        "${input_fasta}" \\
        checkv_out
    """

    stub:
    """
    mkdir -p checkv_out
    touch checkv_out/quality_summary.tsv checkv_out/completeness.tsv checkv_out/contamination.tsv
    """
}
