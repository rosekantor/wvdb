/*
 * modules/cluster_trim_tools.nf
 * Processes for cluster-based trimming (step 3).
 *
 * Process list:
 *   PARSE_CLUSTERS_ALL      — rank123 mode: all pair files + candidate IDs
 *   EXTRACT_TRIMMING_SEQS   — seqkit grep: rank1/2/3 sequences into candidates.fasta
 *   TRIM_GENOMES            — nucmer-based trimming (called 3× via aliases)
 *   PREPARE_CHECKV_INPUT    — cat trim12+trim13 with suffixed IDs for CheckV
 *   PICK_BEST_TRIM1         — split by CheckV completeness after trim1
 *   FETCH_BLAST_INPUT_SEQS  — retrieve untrimmed seqs for 2-member cluster failures
 *   FILTER_PAIRS23          — restrict pairs23.tsv to incomplete rank1 clusters
 *   CHECKV                  — checkv end_to_end (called twice via aliases)
 *   COMPLETENESS_FILTER     — completeness_filter.py (called twice via aliases)
 *
 * publishDir paths use task.ext.publish_dir set via nextflow.config withName selectors.
 *
 * Tools required in PATH:
 *   parse_clusters.py, trim_genomes.py, pick_best_trim1.py,
 *   completeness_filter.py, seqkit, nucmer, show-coords, checkv (via bin/ or conda)
 */

// ---------------------------------------------------------------------------
// PARSE_CLUSTERS_ALL
// Single pass over vclust clusters producing all pair files and candidate IDs.
// ---------------------------------------------------------------------------
process PARSE_CLUSTERS_ALL {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path clusters   // vclust_clusters.tsv
    path lengths    // filtered_all.length.txt

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
// Single seqkit grep extracting all rank1/2/3 sequences.
// The resulting candidates.fasta is pre-indexed by trim_genomes.py (seqkit faidx)
// so all parallel nucmer workers can extract sequences without rebuilding the index.
// ---------------------------------------------------------------------------
process EXTRACT_TRIMMING_SEQS {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path candidate_ids   // candidate_ids.txt from PARSE_CLUSTERS_ALL
    path all_fasta       // filtered_all.fasta

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
// Nucmer-based trimming. Called three times via aliases in cluster_trim.nf:
//   TRIM_GENOMES     → trim12 (rank1 vs rank2, all clusters)
//   TRIM_GENOMES_13  → trim13 (rank1 vs rank3, 3+-member clusters)
//   TRIM_GENOMES_23  → trim23 (rank2 vs rank3, incomplete rank1 subset only)
//
// tag_label scopes output filenames to avoid collisions between parallel calls.
// trim_genomes.py pre-builds a seqkit faidx index on candidates_fasta once,
// then all per-pair nucmer extractions use the small pre-indexed file.
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
// PREPARE_CHECKV_INPUT
// Cats trim12 and trim13 FASTAs with suffixed IDs so CheckV can assess both
// trimmings of each rank1 sequence in a single run.
// Output IDs: <rank1_id>_trim12 and <rank1_id>_trim13
// PICK_BEST_TRIM1 uses completeness from both to choose the winner.
// ---------------------------------------------------------------------------
process PREPARE_CHECKV_INPUT {
    label 'cpu_low'

    input:
    path trim12_fasta   // trim12.trimmed.fasta from TRIM_GENOMES
    path trim13_fasta   // trim13.trimmed.fasta from TRIM_GENOMES_13

    output:
    path "checkv_input.fasta", emit: checkv_input_fasta

    script:
    """
    # Add _trim12 / _trim13 suffixes to FASTA headers so CheckV IDs are unique.
    # awk appends the suffix to the sequence ID (first word after '>'),
    # preserving any description text that follows.
    awk '/^>/{sub(/^>/, ""); id=\$1; \$1=""; printf ">%s_trim12%s\n", id, \$0; next} {print}' \
        "${trim12_fasta}" > checkv_input.fasta

    if [[ -s "${trim13_fasta}" ]]; then
        awk '/^>/{sub(/^>/, ""); id=\$1; \$1=""; printf ">%s_trim13%s\n", id, \$0; next} {print}' \
            "${trim13_fasta}" >> checkv_input.fasta
    fi

    [[ -s checkv_input.fasta ]] || { echo "ERROR: checkv_input.fasta is empty" >&2; exit 1; }
    """

    stub:
    """
    touch checkv_input.fasta
    """
}

// ---------------------------------------------------------------------------
// PICK_BEST_TRIM1
// Selects longer of trim12/trim13 per cluster, splits by CheckV completeness.
//   complete_fasta      → step 5 (recluster)
//   incomplete_ids      → FILTER_PAIRS23 → TRIM_GENOMES_23 (branch B)
//   blast_trim_ids      → FETCH_BLAST_INPUT_SEQS → step 4 (blast_trim)
// ---------------------------------------------------------------------------
process PICK_BEST_TRIM1 {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path trim12_fasta     // trim12.trimmed.fasta from TRIM_GENOMES
    path trim13_fasta     // trim13.trimmed.fasta from TRIM_GENOMES_13
    path checkv_tsv       // quality_summary.tsv from CHECKV
    path pairs12          // pairs12.tsv
    path pairs13          // pairs13.tsv

    output:
    path "complete_reps.fasta",      emit: complete_fasta
    path "incomplete_rank1_ids.txt", emit: incomplete_ids   // → branch B (trim23)
    path "blast_trim_ids.txt",       emit: blast_trim_ids   // → step 4 (blast trim)

    script:
    def trim13_arg  = trim13_fasta.name != 'NO_FILE' ? "--trim13-fasta ${trim13_fasta}" : ""
    def pairs13_arg = pairs13.name     != 'NO_FILE'  ? "--pairs13 ${pairs13}"           : ""
    """
    pick_best_trim1.py \\
        --trim12-fasta    "${trim12_fasta}" \\
        ${trim13_arg} \\
        --checkv-tsv      "${checkv_tsv}" \\
        --pairs12         "${pairs12}" \\
        ${pairs13_arg} \\
        --threshold       ${params.completeness} \\
        --complete-fasta  complete_reps.fasta \\
        --incomplete-ids  incomplete_rank1_ids.txt \\
        --blast-trim-ids  blast_trim_ids.txt
    """

    stub:
    """
    touch complete_reps.fasta incomplete_rank1_ids.txt blast_trim_ids.txt
    """
}

// ---------------------------------------------------------------------------
// FETCH_BLAST_INPUT_SEQS
// Retrieves untrimmed sequences for rank1 IDs of 2-member clusters where
// trim12 was incomplete. These are passed to the blast_trim step.
// ---------------------------------------------------------------------------
process FETCH_BLAST_INPUT_SEQS {
    label 'cpu_low'

    input:
    path ids_file      // blast_trim_ids.txt from PICK_BEST_TRIM1
    path all_fasta     // filtered_all.fasta

    output:
    path "blast_trim_from_trim1.fasta", emit: blast_trim_fasta

    script:
    """
    if [[ -s "${ids_file}" ]]; then
        seqkit grep -f "${ids_file}" "${all_fasta}" -o blast_trim_from_trim1.fasta
    else
        touch blast_trim_from_trim1.fasta
    fi
    """

    stub:
    """
    touch blast_trim_from_trim1.fasta
    """
}

// ---------------------------------------------------------------------------
// FILTER_PAIRS23
// Restricts pairs23.tsv to only clusters whose rank1 ID appears in
// incomplete_rank1_ids.txt (output of PICK_BEST_TRIM1).
// This ensures TRIM_GENOMES_23 (nucmer) only runs on clusters that need it,
// rather than all 3+-member clusters.
// ---------------------------------------------------------------------------
process FILTER_PAIRS23 {
    label 'cpu_low'

    input:
    path pairs23           // pairs23.tsv: col1=rank2 query, col2=rank3 target
    path pairs12           // pairs12.tsv: col1=rank1, col2=rank2 (to map rank1→rank2)
    path incomplete_ids    // incomplete_rank1_ids.txt from PICK_BEST_TRIM1

    output:
    path "pairs23_filtered.tsv", emit: pairs23_filtered

    script:
    """
    # Build a lookup: rank2_id → rank1_id from pairs12
    # Then keep pairs23 rows where the corresponding rank1 is in incomplete_ids
    python3 - << 'PYEOF'
import sys

# Load incomplete rank1 IDs
with open("${incomplete_ids}") as f:
    incomplete = {line.strip() for line in f if line.strip()}

# Load pairs12: rank1 → rank2
rank1_to_rank2 = {}
with open("${pairs12}") as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) >= 2:
            rank1_to_rank2[parts[0]] = parts[1]

# Invert: rank2 → rank1
rank2_to_rank1 = {v: k for k, v in rank1_to_rank2.items()}

# Filter pairs23: keep rows where rank2's corresponding rank1 is incomplete
kept = 0
with open("${pairs23}") as fin, open("pairs23_filtered.tsv", "w") as fout:
    for line in fin:
        parts = line.strip().split('\\t')
        if len(parts) < 2:
            continue
        rank2 = parts[0]
        rank1 = rank2_to_rank1.get(rank2)
        if rank1 and rank1 in incomplete:
            fout.write(line)
            kept += 1

print(f"FILTER_PAIRS23: kept {kept} pairs from incomplete rank1 clusters", file=sys.stderr)
PYEOF
    """

    stub:
    """
    touch pairs23_filtered.tsv
    """
}

// ---------------------------------------------------------------------------
// CHECKV
// Runs checkv end_to_end. Called twice in cluster_trim.nf:
//   CHECKV        → trim1 quality assessment
//   CHECKV_TRIM23 → trim23 quality assessment
// task.ext.publish_dir set per-alias in nextflow.config.
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
// COMPLETENESS_FILTER
// Wraps completeness_filter.py. Called twice in cluster_trim.nf:
//   COMPLETENESS_FILTER_TRIM23 → after trim23 checkv
// emit_incomplete_fasta=true writes untrimmed seqs for blast_trim step.
// task.ext.publish_dir set per-alias in nextflow.config.
// ---------------------------------------------------------------------------
process COMPLETENESS_FILTER {
    label 'cpu_low'

    tag "${quality_summary.simpleName}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path quality_summary
    path trimmed_fasta
    path untrimmed_fasta       // filtered_all.fasta
    val  threshold
    val  emit_incomplete_fasta // true: write untrimmed incomplete seqs

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
