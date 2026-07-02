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
// Produces ani.tsv (and optionally aln.tsv).
// Does NOT produce an ids file — that is written by vclust cluster.
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
    path "vclust_ani.ids.tsv", emit: ids_tsv   // required input for vclust cluster

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

    # vclust cluster requires a TSV of sequence IDs and lengths (id, seq_len, no_parts).
    # vclust align writes this internally to a temp dir and discards it, so we
    # generate it here from the input FASTA using seqkit.
    # vclust cluster --ids requires a TSV with columns: id, seq_len, no_parts
    # Generate from the input FASTA (no_parts is always 1 for standard FASTA).
    printf "id\tseq_len\tno_parts\n" > vclust_ani.ids.tsv
    seqkit fx2tab -nl "${fasta}" | awk 'BEGIN{OFS="\t"}{print \$1, \$2, 1}' >> vclust_ani.ids.tsv
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
