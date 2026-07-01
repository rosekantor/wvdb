#!/usr/bin/env nextflow
/*
 * main.nf  —  wvdb_build: viral genome database build pipeline
 *
 * Steps:
 *   1. Collect genomes    (COLLECT_GENOMES)
 *   2. Cluster            (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER)
 *   3. Cluster trim       (CLUSTER_TRIM subworkflow)
 *                            Branch A: parallel trim12 + trim13 → checkv → pick_best
 *                            Branch B: sequential trim23 on incomplete rank1 subset
 *   4. BLAST trim         (BLAST_TRIM subworkflow)
 *                            singletons + cluster-trim failures → BLAST → trim → checkv
 *   5. Recluster          (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER → GET_CENTROIDS)
 *
 * Usage:
 *   nextflow run main.nf -profile slurm [--fastqdir /path] [--outdir /path]
 *   nextflow run main.nf -stub -profile local   # validate DAG without running tools
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Imports
// vclust processes aliased for step 2 and step 5 invocations.
// ---------------------------------------------------------------------------
include { COLLECT_GENOMES                                   } from './modules/collect'
include { VCLUST_PREFILTER                                  } from './modules/vclust'
include { VCLUST_ALIGN                                      } from './modules/vclust'
include { VCLUST_CLUSTER                                    } from './modules/vclust'
include { VCLUST_PREFILTER as VCLUST_PREFILTER_RECLUST      } from './modules/vclust'
include { VCLUST_ALIGN     as VCLUST_ALIGN_RECLUST          } from './modules/vclust'
include { VCLUST_CLUSTER   as VCLUST_CLUSTER_RECLUST        } from './modules/vclust'
include { GET_CENTROIDS                                     } from './modules/vclust'
include { CLUSTER_TRIM                                      } from './subworkflows/cluster_trim'
include { BLAST_TRIM                                        } from './subworkflows/blast_trim'

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow {

    // Parameter validation
    def errors = []
    if (!params.fastqdir)          errors << "  --fastqdir is required"
    if (!params.outdir)            errors << "  --outdir is required"
    if (!params.checkvdb)          errors << "  --checkvdb is required"
    if (!params.refseq_ev_blastdb) errors << "  --refseq_ev_blastdb is required"
    if (!params.imgvr_blastdb)     errors << "  --imgvr_blastdb is required"
    if (errors) {
        log.error "Parameter errors:\n" + errors.join("\n")
        System.exit(1)
    }

    log.info """
    ============================================
     wvdb_build: viral genome database pipeline
    ============================================
     fastqdir    : ${params.fastqdir}
     outdir      : ${params.outdir}
     ani         : ${params.ani}
     qcov        : ${params.qcov}
     completeness: ${params.completeness}%
     threads     : ${params.threads}
    ============================================
    """.stripIndent()

    def check = !workflow.stubRun
    checkvdb          = file(params.checkvdb,           checkIfExists: check)
    // BLAST dbs passed as strings (not file objects) so Nextflow does not stage
    // them — this keeps all index files (.nhr .nin .nsq etc.) accessible to blastn
    refseq_ev_blastdb = params.refseq_ev_blastdb
    imgvr_blastdb     = params.imgvr_blastdb
    if (!workflow.stubRun) {
        if (!file(params.refseq_ev_blastdb).exists())
            error "refseq_ev_blastdb not found: ${params.refseq_ev_blastdb}"
        if (!file(params.imgvr_blastdb).exists())
            error "imgvr_blastdb not found: ${params.imgvr_blastdb}"
    }
    fastqdir          = file(params.fastqdir,           checkIfExists: check)

    // -----------------------------------------------------------------------
    // Step 1 — Collect and merge input genomes
    // -----------------------------------------------------------------------
    COLLECT_GENOMES(fastqdir)

    // -----------------------------------------------------------------------
    // Step 2 — Initial clustering with vclust
    // -----------------------------------------------------------------------
    VCLUST_PREFILTER(
        COLLECT_GENOMES.out.merged_fasta,
        params.ani
    )

    VCLUST_ALIGN(
        COLLECT_GENOMES.out.merged_fasta,
        VCLUST_PREFILTER.out.prefilter,
        params.ani,
        params.qcov
    )

    VCLUST_CLUSTER(
        VCLUST_ALIGN.out.ani_tsv,
        VCLUST_ALIGN.out.ids_tsv,
        params.ani,
        params.qcov
    )

    // -----------------------------------------------------------------------
    // Step 3 — Cluster-based trimming
    //   minimap2 all-vs-all on rank1/2/3 candidates → pick best alignment
    //   per cluster across pairs12/13/23 → CheckV → split by completeness
    // -----------------------------------------------------------------------
    CLUSTER_TRIM(
        VCLUST_CLUSTER.out.clusters,
        COLLECT_GENOMES.out.lengths,
        COLLECT_GENOMES.out.merged_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Step 4 — BLAST-mode trimming
    //   Receives singletons + incomplete sequences from cluster trim step
    // -----------------------------------------------------------------------
    blast_trim_incomplete = CLUSTER_TRIM.out.blast_trim_fasta

    BLAST_TRIM(
        CLUSTER_TRIM.out.singletons,
        blast_trim_incomplete,
        COLLECT_GENOMES.out.merged_fasta,
        refseq_ev_blastdb,
        imgvr_blastdb,
        checkvdb,
        params.completeness
    )

    // -----------------------------------------------------------------------
    // Step 5 — Recluster all complete representatives
    //   Merges complete outputs from step 3 (branches A + B) and step 4
    // -----------------------------------------------------------------------
    recluster_input = CLUSTER_TRIM.out.complete_fasta
        .mix(BLAST_TRIM.out.complete_fasta)
        .collectFile(name: "recluster_input.fasta")

    VCLUST_PREFILTER_RECLUST(
        recluster_input,
        params.ani
    )

    VCLUST_ALIGN_RECLUST(
        recluster_input,
        VCLUST_PREFILTER_RECLUST.out.prefilter,
        params.ani,
        params.qcov
    )

    VCLUST_CLUSTER_RECLUST(
        VCLUST_ALIGN_RECLUST.out.ani_tsv,
        VCLUST_ALIGN_RECLUST.out.ids_tsv,
        params.ani,
        params.qcov
    )

    GET_CENTROIDS(
        VCLUST_CLUSTER_RECLUST.out.clusters,
        recluster_input
    )
}
