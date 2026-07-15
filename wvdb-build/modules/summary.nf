/*
 * modules/summary.nf
 * PIPELINE_SUMMARY — runs at the end of the workflow, collecting key output
 * FASTAs and TSVs from all steps to produce a sequence count report.
 *
 * Outputs:
 *   pipeline_summary.tsv  — machine-readable counts table
 *   pipeline_summary.md   — human-readable markdown report
 *
 * publishDir: params.outdir (top level, alongside step subdirectories)
 */

process PIPELINE_SUMMARY {
    label 'cpu_low'

    publishDir "${params.outdir}", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path all_fasta                // 1_collect/filtered_all.fasta
    path clusters_tsv             // 2_clustered/vclust_clusters.tsv
    path trimming_candidates      // 3_cluster_trim/trimming_candidates.fasta
    path cluster_trim_complete    // 3_cluster_trim/complete_reps.fasta
    path trim_log                 // 3_cluster_trim/trim_selection.tsv
    path blast_trim_input         // 4_blast_trim/blast_trim_input.fasta
    path blast_trim_complete      // 4_blast_trim/cluster_reps_complete.fasta
    path unvalidated_fasta        // unvalidated/unvalidated_genomes.fasta
    path unvalidated_report       // unvalidated/unvalidated_report.tsv
    path recluster_input          // 5_reclustered/recluster_input.fasta
    path centroids                // 5_reclustered/vclust_centroids.fasta

    output:
    path "pipeline_summary.tsv", emit: summary_tsv
    path "pipeline_summary.md",  emit: summary_md

    script:
    def trim_log_arg     = !trim_log.name.startsWith('NO_FILE') ? "--trim-log ${trim_log}" : ""
    def blast_input_arg  = blast_trim_input.name  != 'NO_FILE' ? "--blast-trim-input  ${blast_trim_input}"  : ""
    def blast_compl_arg  = blast_trim_complete.name != 'NO_FILE' ? "--blast-trim-complete ${blast_trim_complete}" : ""
    def unval_fasta_arg  = unvalidated_fasta.name != 'NO_FILE' ? "--unvalidated-fasta  ${unvalidated_fasta}"  : ""
    def unval_report_arg = unvalidated_report.name != 'NO_FILE' ? "--unvalidated-report ${unvalidated_report}" : ""
    """
    pipeline_summary.py \\
        --all-fasta             "${all_fasta}" \\
        --clusters-tsv          "${clusters_tsv}" \\
        --trimming-candidates   "${trimming_candidates}" \\
        --cluster-trim-complete "${cluster_trim_complete}" \\
        ${blast_input_arg} \\
        ${blast_compl_arg} \\
        ${unval_fasta_arg} \\
        ${unval_report_arg} \\
        ${trim_log_arg} \\
        --recluster-input       "${recluster_input}" \\
        --centroids             "${centroids}" \\
        --out-tsv               pipeline_summary.tsv \\
        --out-md                pipeline_summary.md \\
        --run-date              "\$(date +%Y-%m-%d)" \\
        --ani                   ${params.ani} \\
        --completeness          ${params.completeness}
    """

    stub:
    """
    touch pipeline_summary.tsv pipeline_summary.md
    """
}
