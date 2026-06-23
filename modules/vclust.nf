/*
 * modules/vclust.nf
 * Three vclust processes: prefilter → align → cluster
 * Used for both Step 2 (initial clustering) and Step 6 (reclustering).
 * Aliased at the include level in main.nf for the two invocations.
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
    path "vclust_ani.ids.tsv", emit: ids_tsv

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

    touch vclust_ani.ids.tsv
    """

    stub:
    """
    touch vclust_ani.tsv vclust_ani.aln.tsv vclust_ani.ids.tsv
    """
}

// ---------------------------------------------------------------------------
// VCLUST_CLUSTER
// publishDir is intentionally omitted here — val inputs are not available
// in directive evaluation. Outputs are published via collectFile in main.nf
// (step 6) and via the COLLECT_GENOMES publishDir chain (step 2).
// ---------------------------------------------------------------------------
process VCLUST_CLUSTER {
    label 'cpu_medium'

    tag "${ani_tsv.simpleName}"

    input:
    path ani_tsv
    val  ani
    val  qcov

    output:
    path "vclust_clusters.tsv", emit: clusters
    path "vclust_ani.ids.tsv",  emit: ids

    script:
    """
    vclust cluster \\
        -i "${ani_tsv}" \\
        -o vclust_clusters.tsv \\
        --ids vclust_ani.ids.tsv \\
        --algorithm leiden \\
        --metric ani \\
        --ani ${ani} \\
        --qcov ${qcov} \\
        --out-repr
    """

    stub:
    """
    touch vclust_clusters.tsv vclust_ani.ids.tsv
    """
}

// ---------------------------------------------------------------------------
// GET_CENTROIDS  (Step 6 only)
// ---------------------------------------------------------------------------
process GET_CENTROIDS {
    label 'cpu_low'

    publishDir "${params.outdir}/6_reclustered", mode: 'copy', enabled: !workflow.stubRun

    input:
    path clusters
    path fasta

    output:
    path "vclust_centroids.txt",   emit: centroid_ids
    path "vclust_centroids.fasta", emit: centroid_fasta
    path fasta,                     emit: recluster_input  // publishes recluster_input.fasta

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
