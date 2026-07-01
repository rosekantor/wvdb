/*
 * modules/vclust.nf
 * Three vclust processes: prefilter → align → cluster
 * Used for both Step 2 (initial clustering) and Step 6 (reclustering).
 * Aliased at the include level in main.nf for the two invocations.
 *
 * CLI summary (vclust 1.3.1):
 *   vclust prefilter  -i <fasta> -o <filter>
 *   vclust align      -i <fasta> -o <ani.tsv>   --filter <filter>
 *   vclust cluster    -i <ani.tsv> -o <clusters> --ids <id_map>   (--ids is OUTPUT)
 *
 * Tools required in PATH: vclust v1.3.1+
 */

// ---------------------------------------------------------------------------
// VCLUST_PREFILTER
// ---------------------------------------------------------------------------
process VCLUST_PREFILTER {
    label 'cpu_high'

    tag "${fasta.simpleName}"

    input:
    path fasta
    val  min_ident

    output:
    path "vclust_prefilter.txt", emit: prefilter

    script:
    """
    vclust prefilter \\
        -i "${fasta}" \\
        -o vclust_prefilter.txt \\
        --min-ident ${min_ident} \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch vclust_prefilter.txt
    """
}

// ---------------------------------------------------------------------------
// VCLUST_ALIGN
// Produces ani.tsv, aln.tsv, and vclust_ani.ids.tsv.
// vclust align writes vclust_ani.ids.tsv to the current working directory automatically.
// This file is required as input to vclust cluster --ids.
// ---------------------------------------------------------------------------
process VCLUST_ALIGN {
    label 'cpu_high'

    tag "${fasta.simpleName}"

    input:
    path fasta
    path prefilter
    val  ani
    val  qcov

    output:
    path "vclust_ani.tsv",     emit: ani_tsv
    path "vclust_ani.aln.tsv", emit: aln_tsv
    path "vclust_ani.ids.tsv", emit: ids_tsv   // written by vclust align to CWD

    script:
    """
    vclust align \\
        --filter "${prefilter}" \\
        -i "${fasta}" \\
        --out-ani ${ani} \\
        --out-qcov ${qcov} \\
        -o vclust_ani.tsv \\
        --out-aln vclust_ani.aln.tsv \\
        --threads ${task.cpus}

    [[ -f vclust_ani.ids.tsv ]] || { echo "ERROR: vclust align did not produce vclust_ani.ids.tsv" >&2; exit 1; }
    """

    stub:
    """
    touch vclust_ani.tsv vclust_ani.aln.tsv vclust_ani.ids.tsv
    """
}

// ---------------------------------------------------------------------------
// VCLUST_CLUSTER
// --ids is an OUTPUT path: vclust cluster writes the sequence-ID-to-cluster
// mapping there. It is NOT an input and does not need to exist beforehand.
// ---------------------------------------------------------------------------
process VCLUST_CLUSTER {
    label 'cpu_medium'

    tag "${ani_tsv.simpleName}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path ani_tsv    // output of VCLUST_ALIGN
    path ids_tsv    // sequence ID/length file generated in VCLUST_ALIGN
    val  ani
    val  qcov

    output:
    path "vclust_clusters.tsv", emit: clusters

    script:
    """
    vclust cluster \\
        -i "${ani_tsv}" \\
        -o vclust_clusters.tsv \\
        --ids "${ids_tsv}" \\
        --algorithm leiden \\
        --metric ani \\
        --ani ${ani} \\
        --qcov ${qcov} \\
        --out-repr
    """

    stub:
    """
    touch vclust_clusters.tsv
    """
}

// ---------------------------------------------------------------------------
// GET_CENTROIDS  (Step 5: recluster)
// ---------------------------------------------------------------------------
process GET_CENTROIDS {
    label 'cpu_low'

    publishDir "${params.outdir}/5_reclustered", mode: 'copy', enabled: !workflow.stubRun

    input:
    path clusters
    path fasta

    output:
    path "vclust_centroids.txt",   emit: centroid_ids
    path "vclust_centroids.fasta", emit: centroid_fasta
    path fasta,                    emit: recluster_input  // publishes recluster_input.fasta

    script:
    """
    awk -F'\\t' 'NR>1 {print \$2}' "${clusters}" | LC_ALL=C sort -u > vclust_centroids.txt
    seqkit grep -f vclust_centroids.txt "${fasta}" > vclust_centroids.fasta
    """

    stub:
    """
    touch vclust_centroids.txt vclust_centroids.fasta
    """
}
