/*
 * modules/cluster_trim_tools.nf
 * Processes for cluster-based trimming (step 3).
 *
 * Process list:
 *   PARSE_CLUSTERS_ALL    — rank123 mode: all pair files + candidate IDs
 *   EXTRACT_TRIMMING_SEQS — seqkit grep: rank1/2/3 sequences into candidates.fasta
 *   MINIMAP2_PAIRWISE     — all-vs-all minimap2 on trimming candidates
 *   TRIM_FROM_PAF         — pick best alignment per cluster, write BED, seqkit subseq
 *   CHECKV                — checkv end_to_end
 *   PICK_BEST_TRIM        — split by completeness: complete → recluster, else → blast_trim
 *   FETCH_BLAST_INPUT_SEQS — retrieve untrimmed seqs for blast_trim step
 *
 * Tools required in PATH:
 *   parse_clusters.py, trim_from_paf.py, pick_best_trim.py,
 *   seqkit, minimap2, bgzip, checkv (via bin/ or conda)
 */

// ---------------------------------------------------------------------------
// PARSE_CLUSTERS_ALL
// ---------------------------------------------------------------------------
process PARSE_CLUSTERS_ALL {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path clusters
    path lengths

    output:
    path "pairs12.tsv",        emit: pairs12
    path "pairs13.tsv",        emit: pairs13
    path "pairs23.tsv",        emit: pairs23
    path "singletons.txt",     emit: singletons
    path "candidate_ids.txt",  emit: candidate_ids

    script:
    """
    parse_clusters.py \\
        -c "${clusters}" \\
        -l "${lengths}" \\
        -m rank123 \\
        --out-pairs12       pairs12.tsv \\
        --out-pairs13       pairs13.tsv \\
        --out-pairs23       pairs23.tsv \\
        -s                  singletons.txt \\
        --out-candidate-ids candidate_ids.txt
    """

    stub:
    """
    touch pairs12.tsv pairs13.tsv pairs23.tsv singletons.txt candidate_ids.txt
    """
}

// ---------------------------------------------------------------------------
// EXTRACT_TRIMMING_SEQS
// ---------------------------------------------------------------------------
process EXTRACT_TRIMMING_SEQS {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path candidate_ids
    path all_fasta

    output:
    path "trimming_candidates.fasta", emit: candidates_fasta

    script:
    """
    seqkit grep -f "${candidate_ids}" "${all_fasta}" -o trimming_candidates.fasta

    [[ -s trimming_candidates.fasta ]] || {
        echo "ERROR: trimming_candidates.fasta is empty" >&2
        exit 1
    }
    """

    stub:
    """
    touch trimming_candidates.fasta
    """
}

// ---------------------------------------------------------------------------
// MINIMAP2_PAIRWISE
// All-vs-all minimap2 on trimming candidates (rank1/2/3 sequences only).
// asm5 preset: <5% divergence, appropriate for sequences in the same cluster.
// Full PAF published bgzipped for downstream network analysis.
// ---------------------------------------------------------------------------
process MINIMAP2_PAIRWISE {
    label 'cpu_high'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun,
        saveAs: { fn -> fn.endsWith('.gz') ? fn : null }

    input:
    path candidates_fasta

    output:
    path "alignments.paf",    emit: paf
    path "alignments.paf.gz", emit: paf_gz

    script:
    """
    minimap2 \\
        -x asm5 \\
        -c \\
        --cs \\
        -t ${task.cpus} \\
        "${candidates_fasta}" \\
        "${candidates_fasta}" \\
        > alignments.paf

    bgzip -k alignments.paf
    """

    stub:
    """
    touch alignments.paf alignments.paf.gz
    """
}

// ---------------------------------------------------------------------------
// TRIM_FROM_PAF
// Reads PAF + all three pair files, picks best alignment per cluster by
// alignment block length, writes BED, runs seqkit subseq → trimmed.fasta.
// For rank12/13: rank1 is trimmed. For rank23: rank2 is trimmed.
// ---------------------------------------------------------------------------
process TRIM_FROM_PAF {
    label 'cpu_medium'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path paf
    path pairs12
    path pairs13
    path pairs23
    path candidates_fasta

    output:
    path "trimmed.fasta", emit: trimmed_fasta
    path "trimming.bed",  emit: trimming_bed

    script:
    """
    trim_from_paf.py \\
        --paf     "${paf}" \\
        --pairs12 "${pairs12}" \\
        --pairs13 "${pairs13}" \\
        --pairs23 "${pairs23}" \\
        --fasta   "${candidates_fasta}" \\
        --outdir  . \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch trimmed.fasta trimming.bed
    """
}

// ---------------------------------------------------------------------------
// CHECKV
// ---------------------------------------------------------------------------
process CHECKV {
    label 'cpu_medium'

    tag "${trimmed_fasta.simpleName}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

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
// PICK_BEST_TRIM
// Splits trimmed output by CheckV completeness.
// ---------------------------------------------------------------------------
process PICK_BEST_TRIM {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path trimmed_fasta
    path checkv_tsv

    output:
    path "complete_reps.fasta",  emit: complete_fasta
    path "blast_trim_ids.txt",   emit: blast_trim_ids

    script:
    """
    pick_best_trim.py \\
        --trimmed-fasta  "${trimmed_fasta}" \\
        --checkv-tsv     "${checkv_tsv}" \\
        --threshold      ${params.completeness} \\
        --complete-fasta complete_reps.fasta \\
        --blast-trim-ids blast_trim_ids.txt
    """

    stub:
    """
    touch complete_reps.fasta blast_trim_ids.txt
    """
}

// ---------------------------------------------------------------------------
// FETCH_BLAST_INPUT_SEQS
// Retrieves untrimmed sequences for clusters that were incomplete after
// cluster trimming, to pass to the blast_trim step.
// ---------------------------------------------------------------------------
process FETCH_BLAST_INPUT_SEQS {
    label 'cpu_low'

    input:
    path ids_file
    path all_fasta

    output:
    path "blast_trim_from_cluster.fasta", emit: blast_trim_fasta

    script:
    """
    if [[ -s "${ids_file}" ]]; then
        seqkit grep -f "${ids_file}" "${all_fasta}" -o blast_trim_from_cluster.fasta
    else
        touch blast_trim_from_cluster.fasta
    fi
    """

    stub:
    """
    touch blast_trim_from_cluster.fasta
    """
}

// ---------------------------------------------------------------------------
// COMPLETENESS_FILTER
// Used by the blast_trim subworkflow (imported from here).
// ---------------------------------------------------------------------------
process COMPLETENESS_FILTER {
    label 'cpu_low'

    tag "${quality_summary.simpleName}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path quality_summary
    path trimmed_fasta
    path untrimmed_fasta
    val  threshold
    val  emit_incomplete_fasta

    output:
    path "cluster_reps_complete.fasta",       emit: complete_fasta
    path "cluster_reps_incomplete_or_na.txt", emit: incomplete_ids
    path "cluster_reps_incomplete.fasta", optional: true, emit: incomplete_fasta

    script:
    def inc_flag = emit_incomplete_fasta \
        ? "--incomplete-untrimmed-fasta cluster_reps_incomplete.fasta" : ""
    """
    completeness_filter.py \\
        -q "${quality_summary}" \\
        --trimmed-fasta   "${trimmed_fasta}" \\
        --untrimmed-fasta "${untrimmed_fasta}" \\
        --threshold       ${threshold} \\
        -o                cluster_reps_complete.fasta \\
        --incomplete-ids  cluster_reps_incomplete_or_na.txt \\
        ${inc_flag}
    """

    stub:
    """
    touch cluster_reps_complete.fasta cluster_reps_incomplete_or_na.txt
    ${emit_incomplete_fasta ? 'touch cluster_reps_incomplete.fasta' : ''}
    """
}
