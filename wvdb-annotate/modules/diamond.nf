/*
 * modules/diamond.nf
 * DIAMOND — protein similarity search against a reference database.
 * Uses predicted proteins from geNomad (proteins.faa) as query.
 * Scaffolded: DIAMOND is installed in the conda env but downstream
 * parsing scripts are pending. Output is published for manual inspection.
 *
 * diamond_db is passed as val to prevent index file staging issues.
 */

process DIAMOND {
    label 'cpu_high'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path  proteins_faa    // predicted proteins from GENOMAD
    val   diamond_db      // path to prebuilt .dmnd database (val not path)

    output:
    path "diamond_out.tsv", emit: diamond_tsv

    script:
    """
    diamond blastp \\
        --query   "${proteins_faa}" \\
        --db      "${diamond_db}" \\
        --out     diamond_out.tsv \\
        --outfmt  6 qseqid sseqid pident length mismatch gapopen \\
                    qstart qend sstart send evalue bitscore qlen slen \\
        --evalue  1e-5 \\
        --max-target-seqs 10 \\
        --threads ${task.cpus} \\
        --more-sensitive
    """

    stub:
    """
    touch diamond_out.tsv
    """
}
