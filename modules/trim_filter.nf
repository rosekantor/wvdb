/*
 * modules/trim_filter.nf
 * Processes shared by Branch A and Branch B (and their shared subworkflow).
 *
 *   PARSE_CLUSTERS       — wraps parse_clusters_v2.py
 *   TRIM_GENOMES         — wraps trim_genomes.py
 *   CHECKV               — wraps checkv end_to_end
 *   COMPLETENESS_FILTER  — wraps completeness_filter.py
 *
 * Tools required in PATH:
 *   parse_clusters_v2.py, trim_genomes.py, completeness_filter.py (via bin/)
 *   checkv v1.0.3+, nucmer v4.0.1+
 */

// ---------------------------------------------------------------------------
// PARSE_CLUSTERS
// Produces a pairs TSV (for trim_genomes) and a singletons file.
// mode: "rank12" (branch A) or "rank23" (branch B)
// restrict_ids: optional path to incomplete IDs file; pass [] to skip
// ---------------------------------------------------------------------------
process PARSE_CLUSTERS {
    label 'cpu_low'

    tag "${mode}"

    input:
    path clusters       // vclust_clusters.tsv
    path lengths        // filtered_all.length.txt
    val  mode           // "rank12" or "rank23"
    path restrict_ids   // file of IDs to restrict to, or [] for none

    output:
    path "cluster_pairs.tsv",    emit: pairs
    path "cluster_singletons.txt", emit: singletons

    script:
    def restrict_flag = restrict_ids.name != 'NO_FILE' ? "--restrict-reps ${restrict_ids}" : ""
    """
    parse_clusters_v2.py \\
        -c "${clusters}" \\
        -l "${lengths}" \\
        -o cluster_pairs.tsv \\
        -s cluster_singletons.txt \\
        -m ${mode} \\
        ${restrict_flag}
    """

    stub:
    """
    touch cluster_pairs.tsv cluster_singletons.txt
    """
}

// ---------------------------------------------------------------------------
// TRIM_GENOMES
// Cluster-mode trimming: trim_genomes.py -c pairs -f fasta -o outdir
// ---------------------------------------------------------------------------
process TRIM_GENOMES {
    label 'cpu_high'

    tag "${pairs.simpleName}"

    input:
    path pairs          // cluster_pairs.tsv
    path fasta          // all_genomes.fasta

    output:
    path "trimmed.fasta",  emit: trimmed_fasta
    path "trimming.bed",   emit: trimming_bed   // produced by trim_genomes.py

    script:
    """
    trim_genomes.py \\
        -c "${pairs}" \\
        -f "${fasta}" \\
        -o . \\
        -t ${task.cpus}
    """

    stub:
    """
    touch trimmed.fasta trimming.bed
    """
}

// ---------------------------------------------------------------------------
// CHECKV
// Runs checkv end_to_end on a trimmed FASTA.
// ---------------------------------------------------------------------------
process CHECKV {
    label 'cpu_medium'

    tag "${trimmed_fasta.simpleName}"

    input:
    path trimmed_fasta
    path checkvdb

    output:
    path "checkv_out/quality_summary.tsv", emit: quality_summary
    path "checkv_out/",                    emit: checkv_dir

    script:
    """
    checkv end_to_end \\
        -d "${checkvdb}" \\
        -t ${task.cpus} \\
        --remove_tmp \\
        "${trimmed_fasta}" \\
        checkv_out
    """

    stub:
    """
    mkdir -p checkv_out
    touch checkv_out/quality_summary.tsv
    """
}

// ---------------------------------------------------------------------------
// COMPLETENESS_FILTER
// Splits output into complete and incomplete/NA sets.
// emit_incomplete_fasta: true for branch B (needs the untrimmed FASTA for branch C)
// ---------------------------------------------------------------------------
process COMPLETENESS_FILTER {
    label 'cpu_low'

    tag "${quality_summary.simpleName}"

    input:
    path quality_summary        // checkv_out/quality_summary.tsv
    path trimmed_fasta          // trimmed.fasta from TRIM_GENOMES
    path untrimmed_fasta        // filtered_all.fasta (original full set)
    val  threshold              // completeness threshold (%)
    val  emit_incomplete_fasta  // boolean: whether to emit the incomplete FASTA

    output:
    path "cluster_reps_complete.fasta",           emit: complete_fasta
    path "cluster_reps_incomplete_or_na.txt",     emit: incomplete_ids
    path "cluster_reps_incomplete.fasta", optional: true, emit: incomplete_fasta

    script:
    def incomplete_fasta_flag = emit_incomplete_fasta \
        ? "--incomplete-untrimmed-fasta cluster_reps_incomplete.fasta" \
        : ""
    """
    completeness_filter.py \\
        -q "${quality_summary}" \\
        --trimmed-fasta "${trimmed_fasta}" \\
        --untrimmed-fasta "${untrimmed_fasta}" \\
        --threshold ${threshold} \\
        -o cluster_reps_complete.fasta \\
        --incomplete-ids cluster_reps_incomplete_or_na.txt \\
        ${incomplete_fasta_flag}
    """

    stub:
    """
    touch cluster_reps_complete.fasta cluster_reps_incomplete_or_na.txt
    ${emit_incomplete_fasta ? 'touch cluster_reps_incomplete.fasta' : ''}
    """
}
