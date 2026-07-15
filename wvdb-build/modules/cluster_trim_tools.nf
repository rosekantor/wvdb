/*
 * modules/cluster_trim_tools.nf
 * Processes for cluster-based trimming (step 3).
 *
 * Process list:
 *   PARSE_CLUSTERS_ALL    — rank123 mode: all pair files + candidate IDs
 *   EXTRACT_TRIMMING_SEQS — seqkit grep: rank1/2/3 sequences into candidates.fasta
 *   TRIM_GENOMES          — nucmer-based trimming (called 3× via aliases: trim12, trim13, trim23)
 *   MERGE_TRIMMED         — combine all three trim FASTAs with suffixed IDs for CheckV
 *   CHECKV                — checkv end_to_end (one run covers all three trimmings)
 *   PICK_BEST_TRIM        — select best trimming per rank1 by CheckV completeness + length
 *   FETCH_BLAST_INPUT_SEQS — retrieve untrimmed seqs for blast_trim step
 *   COMPLETENESS_FILTER   — used by blast_trim subworkflow (imported from here)
 *
 * Tools required in PATH:
 *   parse_clusters.py, trim_genomes.py, merge_trimmed.py, pick_best_trim.py,
 *   seqkit, nucmer, show-coords, checkv (via bin/ or conda)
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
// TRIM_GENOMES
// Nucmer-based trimming using pre-extracted trimming_candidates.fasta.
// Called three times via aliases in cluster_trim.nf:
//   TRIM_GENOMES     → trim12 (rank1 vs rank2, all clusters)
//   TRIM_GENOMES_13  → trim13 (rank1 vs rank3, 3+-member clusters)
//   TRIM_GENOMES_23  → trim23 (rank2 vs rank3, incomplete rank1 subset only)
//
// trim_genomes.py pre-builds a seqkit faidx index on candidates_fasta once,
// then per-pair nucmer extractions use the small pre-indexed file.
// tag_label scopes output filenames to avoid collisions between parallel calls.
// ---------------------------------------------------------------------------
process TRIM_GENOMES {
    label 'cpu_high'

    tag "${tag_label}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${tag_label}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path pairs_tsv    // pairs12/13/23.tsv
    path candidates   // trimming_candidates.fasta (pre-indexed by trim_genomes.py)
    val  tag_label    // "trim12", "trim13", or "trim23"

    output:
    path "${tag_label}.trimmed.fasta", emit: trimmed_fasta
    path "${tag_label}.trimming.bed",  emit: trimming_bed

    script:
    """
    trim_genomes.py \\
        -c "${pairs_tsv}" \\
        -f "${candidates}" \\
        -o . \\
        -t ${task.cpus}

    mv trimmed.fasta ${tag_label}.trimmed.fasta
    mv trimming.bed  ${tag_label}.trimming.bed
    """

    stub:
    """
    touch ${tag_label}.trimmed.fasta ${tag_label}.trimming.bed
    """
}

// ---------------------------------------------------------------------------
// MERGE_TRIMMED
// Combines trim12, trim13, trim23 FASTAs into one FASTA with suffixed IDs
// (<rank1_id>_trim12 etc.) so CheckV completeness can be attributed to each
// trimming source in PICK_BEST_TRIM.
// ---------------------------------------------------------------------------
process MERGE_TRIMMED {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path trim12_fasta
    path trim13_fasta
    path trim23_fasta
    path pairs12          // needed to translate trim23 rank2 IDs → rank1

    output:
    path "all_trimmed.fasta", emit: all_trimmed_fasta

    script:
    def trim13_arg = trim13_fasta.name.startsWith('NO_FILE') ? "" : "--trim13 \"${trim13_fasta}\""
    def trim23_arg = trim23_fasta.name.startsWith('NO_FILE') ? "" : "--trim23 \"${trim23_fasta}\""
    """
    merge_trimmed.py \\
        --trim12   "${trim12_fasta}" \\
        ${trim13_arg} \\
        ${trim23_arg} \\
        --pairs12  "${pairs12}" \\
        --out      all_trimmed.fasta
    """

    stub:
    """
    touch all_trimmed.fasta
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
// Selects the best trimming per rank1 cluster from trim12/13/23 results,
// using CheckV completeness as primary criterion and length as tiebreaker.
// Splits output into complete → step 5 and incomplete → step 4 (blast_trim).
// ---------------------------------------------------------------------------
process PICK_BEST_TRIM {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path trim12_fasta
    path trim13_fasta
    path trim23_fasta
    path checkv_tsv
    path pairs12          // needed to translate trim23 rank2 IDs → rank1

    output:
    path "complete_reps.fasta",   emit: complete_fasta
    path "blast_trim_ids.txt",    emit: blast_trim_ids
    path "trim_selection.tsv",    emit: trim_log

    script:
    def trim13_arg = trim13_fasta.name.startsWith('NO_FILE') ? "" : "--trim13 \"${trim13_fasta}\""
    def trim23_arg = trim23_fasta.name.startsWith('NO_FILE') ? "" : "--trim23 \"${trim23_fasta}\""
    """
    pick_best_trim.py \\
        --trim12    "${trim12_fasta}" \\
        ${trim13_arg} \\
        ${trim23_arg} \\
        --checkv    "${checkv_tsv}" \\
        --pairs12   "${pairs12}" \\
        --threshold ${params.completeness} \\
        --complete  complete_reps.fasta \\
        --blast-ids blast_trim_ids.txt \\
        --trim-log  trim_selection.tsv
    """

    stub:
    """
    touch complete_reps.fasta blast_trim_ids.txt trim_selection.tsv
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

