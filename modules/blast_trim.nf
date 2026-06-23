/*
 * modules/blast_trim.nf
 * Processes specific to Branch C: singleton + incomplete genome trimming
 * via BLAST-mode reference databases.
 *
 * Process list:
 *   COLLECT_BRANCHC    — seqkit grep singletons + cat with B incompletes
 *   BLASTN             — blastn search against a reference db
 *   BLASTANI           — compute ANI from blastn tabular output
 *   TRIM_GENOMES_BLAST — blast-mode trim_genomes.py (uses -a ani_tsv -d blastdb)
 *   MERGE_BRANCHC      — comm-based dedup: refseq_ev hits + metavr-only additions
 *
 * publishDir paths use task.ext.publish_dir, set by the calling subworkflow.
 *
 * Tools required in PATH:
 *   seqkit v2.9.0+, blastn v2.16.0+
 *   blastani_nayfach.py, trim_genomes.py (via bin/)
 */

// ---------------------------------------------------------------------------
// COLLECT_BRANCHC
// ---------------------------------------------------------------------------
process COLLECT_BRANCHC {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path singletons_ids       // cluster_singletons.txt from PARSE_CLUSTERS (branch A)
    path incomplete_fasta     // cluster_reps_incomplete.fasta from BRANCH_B
    path all_fasta            // filtered_all.fasta

    output:
    path "branchC_genomes.fasta", emit: branchC_fasta

    script:
    """
    seqkit grep -f "${singletons_ids}" "${all_fasta}" -o singletons.fasta
    cat singletons.fasta "${incomplete_fasta}" > branchC_genomes.fasta

    [[ -s branchC_genomes.fasta ]] || { echo "ERROR: branchC_genomes.fasta is empty" >&2; exit 1; }
    """

    stub:
    """
    touch branchC_genomes.fasta
    """
}

// ---------------------------------------------------------------------------
// BLASTN
// ---------------------------------------------------------------------------
process BLASTN {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path query_fasta          // branchC_genomes.fasta
    path blastdb              // BLAST db (staged directory or .fna + index files)
    val  db_name              // short label: "refseq_ev" or "metavr"

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
// TRIM_GENOMES_BLAST
// Blast-mode trim_genomes.py: uses -a ani_tsv and -d blastdb instead of -c pairs.
// ---------------------------------------------------------------------------
process TRIM_GENOMES_BLAST {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path ani_tsv              // output of BLASTANI
    path query_fasta          // branchC_genomes.fasta
    path blastdb              // same BLAST db used in BLASTN
    val  db_name              // "refseq_ev" or "metavr"

    output:
    path "${db_name}.trimmed.fasta", emit: trimmed_fasta
    path "${db_name}.trimming.bed",  emit: trimming_bed

    script:
    """
    trim_genomes.py \\
        -a "${ani_tsv}" \\
        -f "${query_fasta}" \\
        -d "${blastdb}" \\
        -o . \\
        -t ${task.cpus}

    # Rename to db-scoped names to avoid file collisions in MERGE_BRANCHC
    mv trimmed.fasta ${db_name}.trimmed.fasta
    mv trimming.bed  ${db_name}.trimming.bed
    """

    stub:
    """
    touch ${db_name}.trimmed.fasta ${db_name}.trimming.bed
    """
}

// ---------------------------------------------------------------------------
// MERGE_BRANCHC
// Deduplicates trimmed outputs: all refseq_ev + metavr-only additions.
// ---------------------------------------------------------------------------
process MERGE_BRANCHC {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path refseq_ev_trimmed    // *.trimmed.fasta from TRIM_GENOMES_BLAST (refseq_ev)
    path refseq_ev_bed        // *.trimming.bed  from TRIM_GENOMES_BLAST (refseq_ev)
    path metavr_trimmed       // *.trimmed.fasta from TRIM_GENOMES_BLAST (metavr)
    path metavr_bed           // *.trimming.bed  from TRIM_GENOMES_BLAST (metavr)

    output:
    path "branchC_trimmed.fasta", emit: trimmed_fasta

    script:
    """
    awk '{print \$1}' "${refseq_ev_bed}" | LC_ALL=C sort > refseq_ev_hits_sorted.txt
    awk '{print \$1}' "${metavr_bed}"    | LC_ALL=C sort > metavr_hits_sorted.txt

    # IDs trimmed by metavr but NOT by refseq_ev
    LC_ALL=C comm -13 refseq_ev_hits_sorted.txt metavr_hits_sorted.txt > metavr_only_ids.txt

    seqkit grep -f metavr_only_ids.txt "${metavr_trimmed}" > metavr_only.fasta

    cat "${refseq_ev_trimmed}" metavr_only.fasta > branchC_trimmed.fasta

    [[ -s branchC_trimmed.fasta ]] || { echo "ERROR: branchC_trimmed.fasta is empty" >&2; exit 1; }
    """

    stub:
    """
    touch branchC_trimmed.fasta
    """
}
