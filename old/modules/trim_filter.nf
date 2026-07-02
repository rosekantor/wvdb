/*
 * modules/trim_filter.nf
 * Processes for the cluster trimming pipeline (steps 3/4) and shared
 * by branch C (CHECKV, COMPLETENESS_FILTER).
 *
 * Process list:
 *   PARSE_CLUSTERS_ALL    — rank123 mode: produces all pair files + candidate IDs
 *   EXTRACT_TRIMMING_SEQS — seqkit grep: extracts rank1/2/3 sequences
 *   MINIMAP2_PAIRWISE     — single minimap2 all-vs-all on trimming candidates
 *   TRIM_GENOMES          — PAF-based trim_genomes.py (one alias per pair type)
 *   MERGE_BEST_TRIM1      — cat trim12 + trim13, deduplicate by longest per rank1
 *   PICK_BEST_TRIM1       — pick_best_trim1.py: split by completeness
 *   FILTER_TRIM23         — seqkit grep: restrict trim23 to incomplete rank1 IDs
 *   CHECKV                — checkv end_to_end (called twice: trim1, trim23)
 *   COMPLETENESS_FILTER   — completeness_filter.py (called twice: trim1, trim23)
 *
 * publishDir paths use task.ext.publish_dir set via nextflow.config withName selectors.
 *
 * Tools required in PATH (via bin/ or conda env):
 *   parse_clusters_v2.py, trim_genomes.py, pick_best_trim1.py,
 *   completeness_filter.py, seqkit, minimap2, checkv
 */

// ---------------------------------------------------------------------------
// PARSE_CLUSTERS_ALL
// Single pass producing all pair files and candidate ID list for minimap2.
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
    parse_clusters_v2.py \\
        -c "${clusters}" \\
        -l "${lengths}" \\
        -m rank123 \\
        --out-pairs12      pairs12.tsv \\
        --out-pairs13      pairs13.tsv \\
        --out-pairs23      pairs23.tsv \\
        -s                 singletons.txt \\
        --out-candidate-ids candidate_ids.txt
    """

    stub:
    """
    touch pairs12.tsv pairs13.tsv pairs23.tsv singletons.txt candidate_ids.txt
    """
}

// ---------------------------------------------------------------------------
// EXTRACT_TRIMMING_SEQS
// Single seqkit grep extracting all rank1/2/3 sequences for minimap2.
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
        echo "ERROR: trimming_candidates.fasta is empty — check candidate_ids.txt" >&2
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
// All-vs-all minimap2 on trimming candidates.
// asm5 preset: designed for sequences with <5% divergence (appropriate for
// sequences within a cluster at 95% ANI).
// Full PAF is published bgzipped for downstream network analysis.
// ---------------------------------------------------------------------------
process MINIMAP2_PAIRWISE {
    label 'cpu_high'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun,
        saveAs: { fn -> fn.endsWith('.gz') ? fn : null }  // only publish the bgzipped PAF

    input:
    path candidates_fasta   // trimming_candidates.fasta

    output:
    path "alignments.paf",    emit: paf          // plain PAF for downstream processes
    path "alignments.paf.gz", emit: paf_gz       // bgzipped PAF published for network analysis

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
// TRIM_GENOMES
// PAF-based trimming. Called three times (trim12, trim13, trim23) via aliases.
// tag_label distinguishes the three calls in logs and publishDir.
// ---------------------------------------------------------------------------
process TRIM_GENOMES {
    label 'cpu_medium'

    tag "${tag_label}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${tag_label}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path paf          // alignments.paf from MINIMAP2_PAIRWISE
    path pairs_tsv    // pairs12/13/23.tsv
    path candidates   // trimming_candidates.fasta
    val  tag_label    // "trim12", "trim13", or "trim23"

    output:
    path "${tag_label}.trimmed.fasta", emit: trimmed_fasta
    path "${tag_label}.trimming.bed",  emit: trimming_bed

    script:
    """
    trim_genomes.py \\
        -p "${paf}" \\
        --pairs "${pairs_tsv}" \\
        -f "${candidates}" \\
        -o . \\
        -t ${task.cpus}

    # Rename to tag-scoped names to avoid file collisions in MERGE_BEST_TRIM1
    mv trimmed.fasta ${tag_label}.trimmed.fasta
    mv trimming.bed  ${tag_label}.trimming.bed
    """

    stub:
    """
    touch ${tag_label}.trimmed.fasta ${tag_label}.trimming.bed
    """
}

// ---------------------------------------------------------------------------
// MERGE_BEST_TRIM1
// Merges trim12 and trim13 FASTAs, keeping the longer sequence per rank1 ID.
// The merged output is the input to CHECKV_TRIM1.
// ---------------------------------------------------------------------------
process MERGE_BEST_TRIM1 {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path trim12_fasta   // trimmed.fasta from TRIM_GENOMES (trim12)
    path trim13_fasta   // trimmed.fasta from TRIM_GENOMES (trim13); may be empty

    output:
    path "best_trim1.fasta", emit: best_trim1_fasta

    script:
    """
    # If trim13 is empty (all 2-member clusters), just use trim12 directly
    if [[ ! -s "${trim13_fasta}" ]]; then
        cp "${trim12_fasta}" best_trim1.fasta
    else
        # Concatenate both; pick_best_trim1.py will select the longer per rank1 ID
        cat "${trim12_fasta}" "${trim13_fasta}" > best_trim1_all.fasta
        # Use seqkit to deduplicate keeping the longest per ID
        seqkit rmdup -s best_trim1_all.fasta 2>/dev/null || true
        # Length-based selection: sort by length desc, then keep first occurrence of each ID
        seqkit sort -lr best_trim1_all.fasta | seqkit rmdup -n -o best_trim1.fasta
    fi

    [[ -s best_trim1.fasta ]] || {
        echo "ERROR: best_trim1.fasta is empty" >&2
        exit 1
    }
    """

    stub:
    """
    touch best_trim1.fasta
    """
}

// ---------------------------------------------------------------------------
// PICK_BEST_TRIM1
// Selects longer of trim12/trim13 per cluster, splits by CheckV completeness.
// ---------------------------------------------------------------------------
process PICK_BEST_TRIM1 {
    label 'cpu_low'

    publishDir "${params.outdir}/3_cluster_trim", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path trim12_fasta     // trimmed.fasta from TRIM_GENOMES (trim12)
    path trim13_fasta     // trimmed.fasta from TRIM_GENOMES (trim13); may be empty
    path checkv_tsv       // quality_summary.tsv from CHECKV_TRIM1
    path pairs12          // pairs12.tsv
    path pairs13          // pairs13.tsv; may be empty

    output:
    path "complete_reps.fasta",        emit: complete_fasta
    path "incomplete_rank1_ids.txt",   emit: incomplete_ids     // → try trim23
    path "branch_c_ids.txt",           emit: branch_c_ids       // → branch C directly

    script:
    def trim13_arg = trim13_fasta.name != 'NO_FILE' ? "--trim13-fasta ${trim13_fasta}" : ""
    def pairs13_arg = pairs13.name != 'NO_FILE'     ? "--pairs13 ${pairs13}"           : ""
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
        --branch-c-ids    branch_c_ids.txt
    """

    stub:
    """
    touch complete_reps.fasta incomplete_rank1_ids.txt branch_c_ids.txt
    """
}

// ---------------------------------------------------------------------------
// FILTER_TRIM23
// Restricts trim23 FASTA to sequences whose rank1 ID is in incomplete_ids.
// Avoids running CheckV on trim23 for clusters already resolved at trim1.
// ---------------------------------------------------------------------------
process FILTER_TRIM23 {
    label 'cpu_low'

    input:
    path trim23_fasta      // trimmed.fasta from TRIM_GENOMES (trim23)
    path incomplete_ids    // incomplete_rank1_ids.txt from PICK_BEST_TRIM1
    path pairs23           // pairs23.tsv (to map rank2 IDs → rank1 IDs)

    output:
    path "trim23_filtered.fasta", emit: filtered_fasta

    script:
    """
    # pairs23.tsv columns: rank2 (query/trimmed ID), rank3 (target)
    # incomplete_rank1_ids.txt contains rank1 IDs
    # We need to find rank2 IDs whose corresponding rank1 is in incomplete_ids
    # pairs12.tsv maps rank1→rank2; but we have pairs23 mapping rank2→rank3.
    # Use awk to extract rank2 IDs from pairs23 where rank1 is incomplete.
    # Note: trim_genomes.py preserves query (rank2) ID in the trimmed FASTA.

    # Get all rank2 IDs that correspond to incomplete rank1s via pairs12 lookup.
    # Since pairs12 maps rank1→rank2 and pairs23 maps rank2→rank3,
    # we can identify rank2 by cross-referencing through the cluster structure.
    # Simplest: grep trim23_fasta for sequence IDs that are query in pairs23
    # where the pairs23 rank2 ID maps to an incomplete rank1 in pairs12.

    # Write rank2 IDs whose rank1 is incomplete:
    # pairs12: col1=rank1, col2=rank2
    # Extract rank2 values where rank1 is in incomplete_ids
    awk 'NR==FNR{ids[\$1]=1; next} \$1 in ids {print \$2}' \\
        "${incomplete_ids}" <(cat "${pairs23}" | awk '{print \$1, \$2}' OFS='\\t') \\
        > rank2_to_keep.txt 2>/dev/null || true

    # If no rank2 IDs to keep, write empty output
    if [[ ! -s rank2_to_keep.txt ]]; then
        touch trim23_filtered.fasta
    else
        seqkit grep -f rank2_to_keep.txt "${trim23_fasta}" -o trim23_filtered.fasta
    fi
    """

    stub:
    """
    touch trim23_filtered.fasta
    """
}

// ---------------------------------------------------------------------------
// CHECKV
// Runs checkv end_to_end. Called twice: once for trim1, once for trim23.
// task.ext.publish_dir set via nextflow.config withName selectors.
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
// Wraps completeness_filter.py. Called twice: after trim1 and after trim23.
// emit_incomplete_fasta=true writes untrimmed sequences for branch C.
// task.ext.publish_dir set via nextflow.config withName selectors.
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
    val  emit_incomplete_fasta // true: write untrimmed incomplete seqs for branch C

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

// ---------------------------------------------------------------------------
// FETCH_BRANCH_C_SEQS
// Retrieves untrimmed sequences for rank1 IDs that failed cluster trimming
// (2-member clusters where trim12 was incomplete). These are sent to branch C.
// ---------------------------------------------------------------------------
process FETCH_BRANCH_C_SEQS {
    label 'cpu_low'

    input:
    path ids_file      // branch_c_ids.txt from PICK_BEST_TRIM1
    path all_fasta     // filtered_all.fasta

    output:
    path "branch_c_from_trim1.fasta", emit: branch_c_fasta

    script:
    """
    if [[ -s "${ids_file}" ]]; then
        seqkit grep -f "${ids_file}" "${all_fasta}" -o branch_c_from_trim1.fasta
    else
        touch branch_c_from_trim1.fasta
    fi
    """

    stub:
    """
    touch branch_c_from_trim1.fasta
    """
}
