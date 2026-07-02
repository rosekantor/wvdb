/*
 * modules/rdrpcatch.nf
 * RDRPCATCH — RNA-dependent RNA polymerase detection.
 * Scans nucleotide sequences for RdRP hallmark domains.
 *
 * RdRPCATCH requires Python 3.12 and cannot share the wvdb-annotate
 * conda env (Python 3.11). It runs in a dedicated 'rdrpcatch' conda
 * environment specified via params.rdrpcatch_conda_env in nextflow.config.
 *
 * Setup:
 *   mamba env create -f envs/wvdb_rdrpcatch.yml
 *   conda activate rdrpcatch
 *   rdrpcatch databases --destination-dir /path/to/rdrp_catch_db
 *
 * Binary path controlled by params.rdrpcatch_bin (default: "rdrpcatch").
 */

process RDRPCATCH {
    label 'cpu_high'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path input_fasta
    path rdrpcatch_db

    output:
    path "rdrpcatch_out/",                emit: rdrpcatch_dir
    path "rdrpcatch_out/*.tsv",           emit: results_tsv

    script:
    """
    ${params.rdrpcatch_bin} scan \\
        -i "${input_fasta}" \\
        -o rdrpcatch_out \\
        -db-dir "${rdrpcatch_db}" \\
        --cpus ${task.cpus} \\
        -seq_type nuc
    """

    stub:
    """
    mkdir -p rdrpcatch_out
    touch rdrpcatch_out/rdrpcatch_results.tsv
    """
}
