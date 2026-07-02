/*
 * modules/rnavirhost.nf
 * RNAVIRHOST — host prediction for RNA viruses.
 * Requires CheckV quality_summary.tsv and geNomad virus_summary.tsv
 * to first generate an order.csv input file.
 *
 * run_rnavirhost.py (in bin/) handles both the order.csv generation
 * and the RNAVirHost run. Upload that script to finalize this process.
 */

process RNAVIRHOST {
    label 'cpu_medium'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path input_fasta
    path checkv_quality       // checkv_out/quality_summary.tsv
    path genomad_virus_summary  // genomad virus_summary.tsv

    output:
    path "rnavirhost_out/",         emit: rnavirhost_dir
    path "rnavirhost_out/*.tsv",    emit: results_tsv, optional: true

    script:
    """
    run_rnavirhost.py \\
        --fasta          "${input_fasta}" \\
        --checkv-tsv     "${checkv_quality}" \\
        --genomad-tsv    "${genomad_virus_summary}" \\
        --outdir         rnavirhost_out \\
        --threads        ${task.cpus}
    """

    stub:
    """
    mkdir -p rnavirhost_out
    touch rnavirhost_out/rnavirhost_predictions.tsv
    """
}
