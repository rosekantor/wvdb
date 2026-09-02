/*
 * modules/blast_trim_tools.nf
 * Processes for BLAST-mode reference-guided trimming (step 4).
 *
 * Process list:
 *   COLLECT_BLAST_INPUT    — seqkit grep singletons + cat with cluster-trim incompletes
 *   GET_QUERY_IDS          — extract sequence IDs from blast_trim_input.fasta
 *   BLASTN                 — blastn search against a reference db
 *   BLASTANI               — compute ANI from blastn tabular output
 *   SELECT_BEST_BLAST_HIT  — route queries: initial db preferred; secondary if no initial hit
 *   TRIM_GENOMES_BLAST     — nucmer-based trimming from ANI TSV + reference db FASTA
 *   MERGE_TRIMMED          — cat initial + secondary trimmed FASTAs (queries are disjoint)
 *   COLLECT_UNVALIDATED    — gather no-hit + incomplete seqs → unvalidated output
 *
 * publishDir paths use task.ext.publish_dir set via nextflow.config withName selectors.
 *
 * Tools required in PATH:
 *   seqkit v2.9.0+, blastn v2.16.0+
 *   blastani_nayfach.py, trim_genomes.py, select_best_blast_hit.py,
 *   collect_unvalidated.py (via bin/)
 */

// ---------------------------------------------------------------------------
// COLLECT_BLAST_INPUT
// ---------------------------------------------------------------------------
process COLLECT_BLAST_INPUT {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path singletons_ids
    path incomplete_fasta
    path all_fasta

    output:
    path "blast_trim_input.fasta", emit: blast_trim_fasta

    script:
    """
    seqkit grep -f "${singletons_ids}" "${all_fasta}" -o singletons.fasta
    cat singletons.fasta "${incomplete_fasta}" > blast_trim_input.fasta

    [[ -s blast_trim_input.fasta ]] || {
        echo "ERROR: blast_trim_input.fasta is empty" >&2
        exit 1
    }
    """

    stub:
    """
    touch blast_trim_input.fasta
    """
}

// ---------------------------------------------------------------------------
// GET_QUERY_IDS
// Extracts sequence IDs from blast_trim_input.fasta for SELECT_BEST_BLAST_HIT.
// ---------------------------------------------------------------------------
process GET_QUERY_IDS {
    label 'cpu_low'

    input:
    path blast_trim_fasta

    output:
    path "blast_trim_input_ids.txt", emit: query_ids

    script:
    """
    seqkit seq -n "${blast_trim_fasta}" > blast_trim_input_ids.txt
    """

    stub:
    """
    touch blast_trim_input_ids.txt
    """
}

// ---------------------------------------------------------------------------
// BLASTN
// db_name val used for output naming: "initial" or "secondary"
// blastdb is val (not path) so index files are not staged away from the db dir
// ---------------------------------------------------------------------------
process BLASTN {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path query_fasta
    val  blastdb              // BLAST db path as val — index files must be co-located
    val  db_name              // "initial" or "secondary"

    output:
    path "${db_name}.blastn.tsv", emit: blastn_tsv

    script:
    """
    blastn \\
        -query "${query_fasta}" \\
        -db "${blastdb}" \\
        -max_target_seqs 10 \\
        -perc_identity 90 \\
        -num_threads ${task.cpus} \\
        -out "${db_name}.blastn.tsv" \\
        -outfmt "6 std qlen slen"
    """

    stub:
    """
    touch ${db_name}.blastn.tsv
    """
}

// ---------------------------------------------------------------------------
// BLASTANI
// ---------------------------------------------------------------------------
process BLASTANI {
    label 'cpu_low'

    tag "${blastn_tsv.simpleName}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${blastn_tsv.simpleName}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path blastn_tsv

    output:
    path "${blastn_tsv.baseName}.ani.tsv", emit: ani_tsv

    script:
    """
    blastani_nayfach.py \\
        -i "${blastn_tsv}" \\
        -o "${blastn_tsv.baseName}.ani.tsv"
    """

    stub:
    """
    touch ${blastn_tsv.baseName}.ani.tsv
    """
}

// ---------------------------------------------------------------------------
// SELECT_BEST_BLAST_HIT
// Routes each query to the initial or secondary db hit, requiring a
// "qualifying" hit: tcov >= params.blast_trim_min_tcov AND
// pid >= params.blast_trim_min_pid. Initial db is preferred when both
// have a qualifying hit; secondary used only for queries with no
// qualifying initial hit. Queries with no qualifying hit in either db
// (including queries with a low-coverage/low-identity hit that doesn't
// clear the threshold) → no_hit_ids.txt → unvalidated (no_qualifying_hit).
// secondary_ani may be an empty file when run_secondary_blast=false.
// ---------------------------------------------------------------------------
process SELECT_BEST_BLAST_HIT {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path initial_ani_tsv      // ANI TSV from BLASTANI (initial db)
    path secondary_ani_tsv    // ANI TSV from BLASTANI_SECONDARY (may be empty)
    path query_ids            // blast_trim_input_ids.txt from GET_QUERY_IDS

    output:
    path "initial_hits.tsv",     emit: initial_hits
    path "secondary_hits.tsv",   emit: secondary_hits
    path "no_hit_ids.txt",       emit: no_hit_ids

    script:
    """
    select_best_blast_hit.py \\
        --initial-ani   "${initial_ani_tsv}" \\
        --secondary-ani "${secondary_ani_tsv}" \\
        --initial-out   initial_hits.tsv \\
        --secondary-out secondary_hits.tsv \\
        --no-hit-ids    no_hit_ids.txt \\
        --all-query-ids "${query_ids}" \\
        --min-tcov      ${params.blast_trim_min_tcov} \\
        --min-pid       ${params.blast_trim_min_pid}
    """

    stub:
    """
    touch initial_hits.tsv secondary_hits.tsv no_hit_ids.txt
    """
}

// ---------------------------------------------------------------------------
// TRIM_GENOMES_BLAST
// Nucmer-based trimming from ANI TSV + reference db FASTA.
// Called twice (initial, secondary) via aliases in blast_trim.nf.
// db_name scopes output filenames to avoid collisions in MERGE_TRIMMED.
// ---------------------------------------------------------------------------
process TRIM_GENOMES_BLAST {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path ani_tsv
    path query_fasta
    val  blastdb              // BLAST db FASTA path as val
    val  db_name              // "initial" or "secondary"

    output:
    path "${db_name}.trimmed.fasta",           emit: trimmed_fasta
    path "${db_name}.trimming.bed",            emit: trimming_bed
    path "${db_name}.no_nucmer_alignment.txt", emit: no_alignment_ids

    script:
    """
    trim_genomes.py \\
        -a "${ani_tsv}" \\
        -f "${query_fasta}" \\
        -d "${blastdb}" \\
        -o . \\
        -t ${task.cpus} \\
        --min-identity ${params.blast_trim_min_pid - params.trim_identity_buffer}

    mv trimmed.fasta ${db_name}.trimmed.fasta
    mv trimming.bed  ${db_name}.trimming.bed
    mv no_nucmer_alignment_ids.txt ${db_name}.no_nucmer_alignment.txt
    """

    stub:
    """
    touch ${db_name}.trimmed.fasta ${db_name}.trimming.bed ${db_name}.no_nucmer_alignment.txt
    """
}

// ---------------------------------------------------------------------------
// MERGE_TRIMMED
// Simple cat of initial + secondary trimmed FASTAs.
// Queries are guaranteed disjoint by SELECT_BEST_BLAST_HIT so no dedup needed.
// ---------------------------------------------------------------------------
process MERGE_TRIMMED {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path initial_trimmed      // initial.trimmed.fasta
    path secondary_trimmed    // secondary.trimmed.fasta (may be empty)

    output:
    path "blast_trimmed.fasta", emit: trimmed_fasta

    script:
    """
    cat "${initial_trimmed}" > blast_trimmed.fasta
    if [[ -s "${secondary_trimmed}" ]]; then
        cat "${secondary_trimmed}" >> blast_trimmed.fasta
    fi

    [[ -s blast_trimmed.fasta ]] || {
        echo "ERROR: blast_trimmed.fasta is empty" >&2
        exit 1
    }
    """

    stub:
    """
    touch blast_trimmed.fasta
    """
}

// ---------------------------------------------------------------------------
// COLLECT_UNVALIDATED
// Gathers all sequences that could not be validated into a single FASTA
// with a summary report. Sequences are written in their UNTRIMMED form.
// Sources: no BLAST hit + incomplete after CheckV trimming.
// ---------------------------------------------------------------------------
process COLLECT_UNVALIDATED {
    label 'cpu_low'

    publishDir "${params.outdir}/unvalidated", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path no_hit_ids           // no_hit_ids.txt from SELECT_BEST_BLAST_HIT
    path no_alignment_ids     // merged initial+secondary no_nucmer_alignment.txt
    path incomplete_ids       // incomplete IDs from COMPLETENESS_FILTER
    path all_fasta            // filtered_all.fasta

    output:
    path "unvalidated_genomes.fasta", emit: unvalidated_fasta
    path "unvalidated_report.tsv",    emit: unvalidated_report

    script:
    """
    collect_unvalidated.py \\
        --no-hit-ids       "${no_hit_ids}" \\
        --no-alignment-ids "${no_alignment_ids}" \\
        --incomplete-ids   "${incomplete_ids}" \\
        --all-fasta        "${all_fasta}" \\
        --out-fasta        unvalidated_genomes.fasta \\
        --out-report       unvalidated_report.tsv
    """

    stub:
    """
    touch unvalidated_genomes.fasta unvalidated_report.tsv
    """
}
