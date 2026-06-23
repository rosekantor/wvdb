#!/usr/bin/env nextflow
/*
 * main.nf  —  wf3 viral genome clustering pipeline
 *
 * Steps:
 *   1. Collect genomes          (COLLECT_GENOMES)
 *   2. Cluster with vclust      (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER)
 *   3. Branch A trimming        (TRIM_AND_FILTER subworkflow, mode=rank12)
 *   4. Branch B trimming        (TRIM_AND_FILTER subworkflow, mode=rank23, restricted to A incompletes)
 *   5. Branch C blast trimming  (subworkflows/branchC.nf — stubbed, see TODO)
 *   6. Recluster                (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER → GET_CENTROIDS)
 *
 * Usage:
 *   nextflow run main.nf -profile slurm [--fastqdir /path] [--outdir /path]
 *   nextflow run main.nf -stub          # validate DAG without running tools
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Imports
// Step 2 and Step 6 both use the three vclust processes — DSL2 requires
// aliased includes for any process invoked more than once in the same workflow.
// ---------------------------------------------------------------------------
include { COLLECT_GENOMES                                   } from './modules/collect'
include { VCLUST_PREFILTER                                  } from './modules/vclust'
include { VCLUST_ALIGN                                      } from './modules/vclust'
include { VCLUST_CLUSTER                                    } from './modules/vclust'
include { VCLUST_PREFILTER as VCLUST_PREFILTER_RECLUST      } from './modules/vclust'
include { VCLUST_ALIGN     as VCLUST_ALIGN_RECLUST          } from './modules/vclust'
include { VCLUST_CLUSTER   as VCLUST_CLUSTER_RECLUST        } from './modules/vclust'
include { GET_CENTROIDS                                     } from './modules/vclust'
include { TRIM_AND_FILTER  as BRANCH_A                      } from './subworkflows/trim_and_filter'
include { TRIM_AND_FILTER  as BRANCH_B                      } from './subworkflows/trim_and_filter'
include { BRANCH_C                                         } from './subworkflows/branchC'

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow {

    // Parameter validation — inline here since def blocks cannot appear
    // between includes and workflow.onComplete in DSL2.
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

    // Resolve database paths so Nextflow can stage them on remote executors
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
        params.ani,
        params.qcov
    )

    // -----------------------------------------------------------------------
    // Step 3 — Branch A: trim rank-1 vs rank-2 contig pairs
    // -----------------------------------------------------------------------
    BRANCH_A(
        VCLUST_CLUSTER.out.clusters,
        COLLECT_GENOMES.out.lengths,
        COLLECT_GENOMES.out.merged_fasta,
        checkvdb,
        "rank12",
        file('NO_FILE'),        // no --restrict-reps for branch A
        params.completeness,
        false                   // don't emit incomplete FASTA
    )

    // -----------------------------------------------------------------------
    // Step 4 — Branch B: trim rank-2 vs rank-3, restricted to A incompletes
    // -----------------------------------------------------------------------
    BRANCH_B(
        VCLUST_CLUSTER.out.clusters,
        COLLECT_GENOMES.out.lengths,
        COLLECT_GENOMES.out.merged_fasta,
        checkvdb,
        "rank23",
        BRANCH_A.out.incomplete_ids,   // --restrict-reps from branch A output
        params.completeness,
        true                           // emit incomplete FASTA for branch C
    )

    // -----------------------------------------------------------------------
    // Step 5 — Branch C: singletons + B incompletes → BLAST-mode trimming
    // -----------------------------------------------------------------------
    BRANCH_C(
        BRANCH_A.out.singletons,
        BRANCH_B.out.incomplete_fasta,
        COLLECT_GENOMES.out.merged_fasta,
        refseq_ev_blastdb,
        imgvr_blastdb,
        checkvdb,
        params.completeness
    )
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // Step 6 — Recluster: merge all complete sets, run vclust again
    // Note: the original script has a manual checkpoint here. In Nextflow this
    // is modelled as a straight cat + recluster. Add an external review step
    // outside the pipeline if manual inspection is still needed.
    // -----------------------------------------------------------------------
    // collectFile merges branch outputs into a single FASTA for reclustering.
    // recluster_input.fasta is published via GET_CENTROIDS publishDir in step 6.
    recluster_input = BRANCH_A.out.complete_fasta
        .mix(BRANCH_B.out.complete_fasta)
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
        params.ani,
        params.qcov
    )

    GET_CENTROIDS(
        VCLUST_CLUSTER_RECLUST.out.clusters,
        recluster_input
    )
}
