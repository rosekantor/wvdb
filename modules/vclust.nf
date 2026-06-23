/*
 * modules/vclust.nf
 * Three vclust processes: prefilter → align → cluster
 * Used for both Step 2 (initial clustering) and Step 6 (reclustering).
 * Parameterised via inputs so the same process definitions serve both steps.
 *
 * Tools required in PATH: vclust v1.3.1+
 */

// ---------------------------------------------------------------------------
// VCLUST_PREFILTER
// Generates a sparse candidate-pair list, filtering on minimum identity.
// ---------------------------------------------------------------------------
process VCLUST_PREFILTER {
    label 'cpu_high'

    tag "${fasta.simpleName}"

    input:
    path fasta          // input FASTA (all genomes or recluster input)
    val  min_ident      // passed as params.ani from the calling workflow

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
// Pairwise alignment of candidate pairs from prefilter, producing ANI table.
// ---------------------------------------------------------------------------
process VCLUST_ALIGN {
    label 'cpu_high'

    tag "${fasta.simpleName}"

    input:
    path fasta          // same FASTA used in prefilter
    path prefilter      // output of VCLUST_PREFILTER
    val  ani            // ANI threshold (e.g. 0.95)
    val  qcov           // query coverage threshold (e.g. 0.85)

    output:
    path "vclust_ani.tsv",     emit: ani_tsv
    path "vclust_ani.aln.tsv", emit: aln_tsv
    path "vclust_ani.ids.tsv", emit: ids_tsv   // written by --ids; used downstream

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

    # vclust align does not write ids itself; create placeholder so VCLUST_CLUSTER
    # can receive it as a declared output even if empty.
    touch vclust_ani.ids.tsv
    """

    stub:
    """
    touch vclust_ani.tsv vclust_ani.aln.tsv vclust_ani.ids.tsv
    """
}

// ---------------------------------------------------------------------------
// VCLUST_CLUSTER
// Leiden clustering from ANI table → cluster assignments + representative IDs.
// ---------------------------------------------------------------------------
process VCLUST_CLUSTER {
    label 'cpu_medium'

    tag "${ani_tsv.simpleName}"

    publishDir "${params.outdir}/${step_dir}", mode: 'copy', pattern: 'vclust_clusters.tsv'
    publishDir "${params.outdir}/${step_dir}", mode: 'copy', pattern: 'vclust_ani.ids.tsv'

    input:
    path ani_tsv        // output of VCLUST_ALIGN
    val  ani            // ANI threshold
    val  qcov           // qcov threshold
    val  step_dir       // subdirectory label, e.g. "2_clustered" or "6_reclustered"

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
// Extract centroid sequences from recluster input FASTA.
// ---------------------------------------------------------------------------
process GET_CENTROIDS {
    label 'cpu_low'

    publishDir "${params.outdir}/6_reclustered", mode: 'copy'

    input:
    path clusters       // vclust_clusters.tsv from VCLUST_CLUSTER (recluster)
    path fasta          // recluster_input.fasta

    output:
    path "vclust_centroids.txt",  emit: centroid_ids
    path "vclust_centroids.fasta",emit: centroid_fasta

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
