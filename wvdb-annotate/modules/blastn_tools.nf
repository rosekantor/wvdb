/*
 * modules/blastn_tools.nf
 * Processes for BLASTn annotation searches.
 *
 * Process list:
 *   BLASTN     — nucleotide search against one database
 *   BLASTANI   — compute pairwise ANI from blastn tabular output
 *
 * Called once per database via a channel of (db_name, db_path) pairs
 * derived from the blastn_dbs CSV in main.nf.
 *
 * blastdb is passed as val (not path) to prevent staging and preserve
 * all index file co-location.
 */

// ---------------------------------------------------------------------------
// BLASTN
// ---------------------------------------------------------------------------
process BLASTN {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path  query_fasta
    val   db_name
    val   db_path         // val not path — index files must remain co-located

    output:
    tuple val(db_name), path("${db_name}.blastn.tsv"), emit: blastn_result

    script:
    """
    blastn \\
        -query    "${query_fasta}" \\
        -db       "${db_path}" \\
        -evalue   ${params.blastn_evalue} \\
        -max_target_seqs ${params.blastn_max_targets} \\
        -num_threads ${task.cpus} \\
        -out      "${db_name}.blastn.tsv" \\
        -outfmt   "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    """

    stub:
    """
    touch ${db_name}.blastn.tsv
    """
}

// ---------------------------------------------------------------------------
// BLASTANI
// Computes pairwise ANI from blastn tabular output.
// Reuses blastani_nayfach.py from bin/.
// ---------------------------------------------------------------------------
process BLASTANI {
    label 'cpu_low'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    tuple val(db_name), path(blastn_tsv)

    output:
    tuple val(db_name), path("${db_name}.blastn.ani.tsv"), emit: ani_result

    script:
    """
    blastani_nayfach.py \\
        -i "${blastn_tsv}" \\
        -o "${db_name}.blastn.ani.tsv"
    """

    stub:
    """
    touch ${db_name}.blastn.ani.tsv
    """
}
