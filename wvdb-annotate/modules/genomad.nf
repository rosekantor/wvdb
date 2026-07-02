/*
 * modules/genomad.nf
 * GENOMAD — virus classification, hallmark gene detection, and gene prediction.
 * Produces predicted proteins used downstream by DIAMOND.
 *
 * --cleanup removes intermediate files to save disk space.
 */

process GENOMAD {
    label 'cpu_high'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path input_fasta
    path genomad_db

    output:
    path "genomad_out/",                                                       emit: genomad_dir
    path "genomad_out/${input_fasta.baseName}_summary/${input_fasta.baseName}_virus_summary.tsv", \
                                                                               emit: virus_summary
    path "genomad_out/${input_fasta.baseName}_find_proviruses/${input_fasta.baseName}_provirus.fna", \
                                                                               emit: provirus_fna,    optional: true
    path "genomad_out/${input_fasta.baseName}_annotate/${input_fasta.baseName}_proteins.faa", \
                                                                               emit: proteins_faa

    script:
    """
    genomad end-to-end \\
        --cleanup \\
        -t ${task.cpus} \\
        "${input_fasta}" \\
        genomad_out \\
        "${genomad_db}"
    """

    stub:
    """
    mkdir -p genomad_out/${input_fasta.baseName}_summary
    mkdir -p genomad_out/${input_fasta.baseName}_find_proviruses
    mkdir -p genomad_out/${input_fasta.baseName}_annotate
    touch genomad_out/${input_fasta.baseName}_summary/${input_fasta.baseName}_virus_summary.tsv
    touch genomad_out/${input_fasta.baseName}_find_proviruses/${input_fasta.baseName}_provirus.fna
    touch genomad_out/${input_fasta.baseName}_annotate/${input_fasta.baseName}_proteins.faa
    """
}
