#!/usr/bin/env nextflow
/*
 * main.nf  —  wf3 viral genome clustering pipeline
 *
 * Steps:
 *   1. Collect genomes          (COLLECT_GENOMES)
 *   2. Cluster with vclust      (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER)
 *   3. Branch A trimming        (TRIM_AND_FILTER subworkflow, mode=rank12)
 *   4. Branch B trimming        (TRIM_AND_FILTER subworkflow, mode=rank23, restricted to A incompletes)
 *   5. Branch C blast trimming  (modules/branchC.nf  — stubbed, see TODO)
 *   6. Recluster                (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER → GET_CENTROIDS)
 *
 * Usage:
 *   nextflow run main.nf -profile slurm [--fastqdir /path] [--outdir /path]
 *   nextflow run main.nf -stub          # validate DAG without running tools
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Imports
// ---------------------------------------------------------------------------
include { COLLECT_GENOMES                        } from './modules/collect'
include { VCLUST_PREFILTER                       } from './modules/vclust'
include { VCLUST_ALIGN                           } from './modules/vclust'
include { VCLUST_CLUSTER                         } from './modules/vclust'
include { GET_CENTROIDS                          } from './modules/vclust'
include { TRIM_AND_FILTER as BRANCH_A            } from './subworkflows/trim_and_filter'
include { TRIM_AND_FILTER as BRANCH_B            } from './subworkflows/trim_and_filter'
// include { BRANCH_C                            } from './subworkflows/branchC'   // TODO: implement in next step

// ---------------------------------------------------------------------------
// Parameter validation (fail fast before any jobs are submitted)
// ---------------------------------------------------------------------------
def validate_params() {
    def errors = []
    if (!params.fastqdir)         errors << "  --fastqdir is required"
    if (!params.outdir)           errors << "  --outdir is required"
    if (!params.checkvdb)         errors << "  --checkvdb is required"
    if (!params.refseq_ev_blastdb) errors << "  --refseq_ev_blastdb is required"
    if (!params.imgvr_blastdb)    errors << "  --imgvr_blastdb is required"
    if (errors) {
        log.error "Parameter errors:\n" + errors.join("\n")
        System.exit(1)
    }
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow {

    validate_params()

    // Print run summary to log
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

    // Resolve database paths as Nextflow path objects so they are staged
    // correctly when running on remote executors.
    checkvdb         = file(params.checkvdb,          checkIfExists: true)
    refseq_ev_blastdb = file(params.refseq_ev_blastdb, checkIfExists: true)
    imgvr_blastdb    = file(params.imgvr_blastdb,     checkIfExists: true)
    fastqdir         = file(params.fastqdir,          checkIfExists: true)

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
        params.qcov,
        "2_clustered"       // publishDir subdirectory label
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
        file('NO_FILE'),       // no restriction on branch A
        params.completeness,
        false                  // don't emit incomplete FASTA for branch A
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
        BRANCH_A.out.incomplete_ids,  // --restrict-reps from branch A
        params.completeness,
        true                   // emit incomplete FASTA so branch C can use it
    )

    // -----------------------------------------------------------------------
    // Step 5 — Branch C: singletons + B incompletes → BLAST-mode trimming
    //
    // TODO: implement subworkflows/branchC.nf then uncomment below.
    // BRANCH_C(
    //     BRANCH_A.out.singletons,
    //     BRANCH_B.out.incomplete_fasta,
    //     COLLECT_GENOMES.out.merged_fasta,
    //     refseq_ev_blastdb,
    //     imgvr_blastdb,
    //     checkvdb,
    //     params.completeness
    // )
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // Step 6 — Recluster: merge all complete sets and recluster
    //
    // Note: the original script has a manual checkpoint comment here.
    // In Nextflow this is modelled as a simple cat + re-run of vclust.
    // Add a manual review step outside the pipeline if needed.
    // -----------------------------------------------------------------------

    // Merge the complete FASTAs from A and B (C will be added once implemented)
    recluster_input = BRANCH_A.out.complete_fasta
        .mix(BRANCH_B.out.complete_fasta)
        // .mix(BRANCH_C.out.complete_fasta)   // uncomment when branchC is ready
        .collectFile(name: "recluster_input.fasta",
                     storeDir: "${params.outdir}/6_reclustered")

    VCLUST_PREFILTER(
        recluster_input,
        params.ani
    )

    VCLUST_ALIGN(
        recluster_input,
        VCLUST_PREFILTER.out.prefilter,
        params.ani,
        params.qcov
    )

    VCLUST_CLUSTER(
        VCLUST_ALIGN.out.ani_tsv,
        params.ani,
        params.qcov,
        "6_reclustered"
    )

    GET_CENTROIDS(
        VCLUST_CLUSTER.out.clusters,
        recluster_input
    )

    log.info "Pipeline complete. Final centroids: ${params.outdir}/6_reclustered/vclust_centroids.fasta"
}

// ---------------------------------------------------------------------------
// Workflow completion handler
// ---------------------------------------------------------------------------
workflow.onComplete {
    log.info (workflow.success
        ? "\nDone! Results in: ${params.outdir}"
        : "\nFailed — check .nextflow.log for details")
}
