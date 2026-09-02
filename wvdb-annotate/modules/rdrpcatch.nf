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
 *
 * RdRPCATCH enforces a hard 300,000 bp sequence length limit and crashes
 * (non-zero exit) on any longer sequence. RNA virus genomes are always
 * far smaller than this, so sequences over the limit are pre-filtered out
 * before the scan and logged to excluded_long_seqs.txt rather than
 * causing the whole process — and pipeline — to fail.
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
    path "excluded_long_seqs.txt",        emit: excluded_ids, optional: true

    script:
    """
    # RdRPCATCH hard-fails on sequences > 300,000 bp. Filter them out first
    # so one oversized vOTU cannot crash the whole scan.
    seqkit seq -M ${params.rdrpcatch_max_seq_len} "${input_fasta}" \\
        > rdrpcatch_input.fasta

    seqkit seq -n -m \$((${params.rdrpcatch_max_seq_len} + 1)) "${input_fasta}" \\
        > excluded_long_seqs.txt || true

    n_excluded=\$(wc -l < excluded_long_seqs.txt)
    if [[ "\$n_excluded" -gt 0 ]]; then
        echo "Warning: \$n_excluded sequence(s) exceed ${params.rdrpcatch_max_seq_len} bp " \\
             "and were excluded from RdRPCATCH scan (see excluded_long_seqs.txt)" >&2
    fi

    ${params.rdrpcatch_bin} scan \\
        -i rdrpcatch_input.fasta \\
        -o rdrpcatch_out \\
        -db-dir "${rdrpcatch_db}" \\
        --cpus ${task.cpus} \\
        -seq_type nuc
    """

    stub:
    """
    mkdir -p rdrpcatch_out
    touch rdrpcatch_out/rdrpcatch_results.tsv excluded_long_seqs.txt
    """
}
