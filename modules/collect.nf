/*
 * modules/collect.nf
 * Step 1: Collect filtered virus and provirus FASTAs from all sample assembly dirs,
 * merge into a single filtered_all.fasta, and compute per-sequence lengths.
 *
 * Expected input directory layout (produced by wf2):
 *   <fastqdir>/<sample>/assembly/<run>/filtered_virus.fasta
 *   <fastqdir>/<sample>/assembly/<run>/filtered_provirus.fasta
 */

process COLLECT_GENOMES {
    label 'cpu_low'

    publishDir "${params.outdir}/1_filtered", mode: 'copy', enabled: !workflow.stubRun

    input:
    path fastqdir   // top-level fastq directory (passed as a path so Nextflow stages it)

    output:
    path "filtered_all.fasta",         emit: merged_fasta
    path "filtered_all.length.txt",    emit: lengths
    path "filtered_virus_all.fasta",   emit: virus_fasta
    path "filtered_provirus_all.fasta",emit: provirus_fasta

    script:
    """
    # Glob for all virus and provirus FASTAs under fastqdir
    virus_files=( "${fastqdir}"/*/assembly/*/filtered_virus.fasta )
    provirus_files=( "${fastqdir}"/*/assembly/*/filtered_provirus.fasta )

    [[ \${#virus_files[@]} -gt 0 ]]    || { echo "ERROR: no filtered_virus.fasta found"    >&2; exit 1; }
    [[ \${#provirus_files[@]} -gt 0 ]] || { echo "ERROR: no filtered_provirus.fasta found" >&2; exit 1; }

    cat "\${virus_files[@]}"    > filtered_virus_all.fasta
    cat "\${provirus_files[@]}" > filtered_provirus_all.fasta
    cat filtered_virus_all.fasta filtered_provirus_all.fasta > filtered_all.fasta

    seqkit fx2tab -nl filtered_all.fasta > filtered_all.length.txt

    [[ -s filtered_all.fasta ]]       || { echo "ERROR: merged FASTA is empty"   >&2; exit 1; }
    [[ -s filtered_all.length.txt ]]  || { echo "ERROR: length file is empty"    >&2; exit 1; }
    """

    stub:
    """
    touch filtered_all.fasta filtered_all.length.txt filtered_virus_all.fasta filtered_provirus_all.fasta
    """
}
