#!/usr/bin/env nextflow
/*
 * main.nf  —  wf3 viral genome clustering pipeline
 *
 * Steps:
 *   1. Collect genomes      (COLLECT_GENOMES)
 *   2. Cluster with vclust  (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER)
 *   3. Cluster trimming     (CLUSTER_TRIM subworkflow)
 *                             parse_clusters → seqkit grep → minimap2 →
 *                             trim12/13/23 → checkv → pick_best → trim23 fallback
 *   4. Branch C             (BRANCH_C subworkflow — BLAST-mode trimming)
 *   5. Recluster            (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER → GET_CENTROIDS)
 *
 * Usage:
 *   nextflow run main.nf -profile slurm [--fastqdir /path] [--outdir /path]
 *   nextflow run main.nf -stub -profile local   # validate DAG without running tools
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Imports
// vclust processes are aliased for the two invocations (step 2 and step 5).
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
include { BRANCH_C                                          } from './subworkflows/branchC'

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
     wf3 viral clustering pipeline
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
    refseq_ev_blastdb = file(params.refseq_ev_blastdb,  checkIfExists: check)
    imgvr_blastdb     = file(params.imgvr_blastdb,      checkIfExists: check)
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
    // Step 3 — Cluster-based trimming (minimap2 pipeline)
    // Replaces separate BRANCH_A and BRANCH_B subworkflows.
    // -----------------------------------------------------------------------
    CLUSTER_TRIM(
        VCLUST_CLUSTER.out.clusters,
        COLLECT_GENOMES.out.lengths,
        COLLECT_GENOMES.out.merged_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Step 4 — Branch C: BLAST-mode trimming
    // Receives:
    //   singletons          — clusters with 1 member (no trimming possible)
    //   branch_c_fasta      — untrimmed rank1 seqs from 2-member clusters
    //                         where trim12 was incomplete
    //   branch_c_fasta23    — untrimmed seqs from clusters where trim23
    //                         was also incomplete (from COMPLETENESS_FILTER_TRIM23)
    // COLLECT_BRANCHC cats these two incomplete FASTAs together.
    // -----------------------------------------------------------------------

    // Merge the two sources of incomplete sequences into one FASTA for branch C
    branch_c_incomplete = CLUSTER_TRIM.out.branch_c_fasta
        .mix(CLUSTER_TRIM.out.branch_c_fasta23)
        .collectFile(name: "branch_c_incomplete.fasta")

    BRANCH_C(
        CLUSTER_TRIM.out.singletons,
        branch_c_incomplete,
        COLLECT_GENOMES.out.merged_fasta,
        refseq_ev_blastdb,
        imgvr_blastdb,
        checkvdb,
        params.completeness
    )

    // -----------------------------------------------------------------------
    // Step 5 — Recluster all complete representatives
    // -----------------------------------------------------------------------
    recluster_input = CLUSTER_TRIM.out.complete_fasta
        .mix(CLUSTER_TRIM.out.complete_fasta23)
        .mix(BRANCH_C.out.complete_fasta)
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
