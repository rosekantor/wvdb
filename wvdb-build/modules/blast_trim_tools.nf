/*
 * modules/blast_trim_tools.nf
 * Processes for BLAST-mode reference-guided trimming (step 4).
 *
 * Process list:
 *   COLLECT_BLAST_INPUT  — seqkit grep singletons + cat with cluster-trim incompletes
 *   BLASTN               — blastn search against a reference db
 *   BLASTANI             — compute ANI from blastn tabular output
 *   TRIM_GENOMES_BLAST   — blast-mode trim_genomes.py (uses -a ani_tsv -d blastdb)
 *   MERGE_BLAST_TRIM     — comm-based dedup: refseq_ev hits + metavr-only additions
 *
 * publishDir paths use task.ext.publish_dir set via nextflow.config withName selectors.
 *
 * Tools required in PATH:
 *   seqkit v2.9.0+, blastn v2.16.0+
 *   blastani_nayfach.py, trim_genomes.py (via bin/)
 */

// ---------------------------------------------------------------------------
// COLLECT_BLAST_INPUT
// Greps singletons from the full FASTA and concatenates with cluster-trim
// incomplete sequences to form the full blast-trim input set.
// ---------------------------------------------------------------------------
process COLLECT_BLAST_INPUT {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path singletons_ids       // cluster_singletons.txt from PARSE_CLUSTERS_ALL
    path incomplete_fasta     // merged incomplete seqs from cluster_trim step
    path all_fasta            // filtered_all.fasta

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
// BLASTN
// ---------------------------------------------------------------------------
process BLASTN {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path query_fasta          // blast_trim_input.fasta
    val  blastdb              // BLAST db path as val — prevents staging so all index files remain accessible
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
// Blast-mode trim_genomes.py: uses -a ani_tsv and -d blastdb.
// db_name scopes output filenames to avoid collisions in MERGE_BLAST_TRIM.
// ---------------------------------------------------------------------------
process TRIM_GENOMES_BLAST {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path ani_tsv              // output of BLASTANI
    path query_fasta          // blast_trim_input.fasta
    val  blastdb              // BLAST db path as val — prevents staging so all index files remain accessible
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

    mv trimmed.fasta ${db_name}.trimmed.fasta
    mv trimming.bed  ${db_name}.trimming.bed
    """

    stub:
    """
    touch ${db_name}.trimmed.fasta ${db_name}.trimming.bed
    """
}

// ---------------------------------------------------------------------------
// MERGE_BLAST_TRIM
// Deduplicates trimmed outputs: keep all refseq_ev trimmed sequences plus
// any metavr-trimmed sequences not already covered by refseq_ev.
// ---------------------------------------------------------------------------
process MERGE_BLAST_TRIM {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path refseq_ev_trimmed    // *.trimmed.fasta from TRIM_GENOMES_BLAST (refseq_ev)
    path refseq_ev_bed        // *.trimming.bed  from TRIM_GENOMES_BLAST (refseq_ev)
    path metavr_trimmed       // *.trimmed.fasta from TRIM_GENOMES_BLAST (metavr)
    path metavr_bed           // *.trimming.bed  from TRIM_GENOMES_BLAST (metavr)

    output:
    path "blast_trim_merged.fasta", emit: trimmed_fasta

    script:
    """
    awk '{print \$1}' "${refseq_ev_bed}" | LC_ALL=C sort > refseq_ev_hits_sorted.txt
    awk '{print \$1}' "${metavr_bed}"    | LC_ALL=C sort > metavr_hits_sorted.txt

    # IDs trimmed by metavr but NOT by refseq_ev
    LC_ALL=C comm -13 refseq_ev_hits_sorted.txt metavr_hits_sorted.txt > metavr_only_ids.txt

    seqkit grep -f metavr_only_ids.txt "${metavr_trimmed}" > metavr_only.fasta

    cat "${refseq_ev_trimmed}" metavr_only.fasta > blast_trim_merged.fasta

    [[ -s blast_trim_merged.fasta ]] || {
        echo "ERROR: blast_trim_merged.fasta is empty" >&2
        exit 1
    }
    """

    stub:
    """
    touch blast_trim_merged.fasta
    """
}
