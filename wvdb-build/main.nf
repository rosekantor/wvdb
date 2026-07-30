#!/usr/bin/env nextflow
/*
 * main.nf  —  wvdb_build: viral genome database build pipeline
 *
 * Steps:
 *   1. Collect genomes    (COLLECT_GENOMES)
 *   2. Cluster            (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER)
 *   3. Cluster trim       (CLUSTER_TRIM subworkflow)
 *                            trim12 ∥ trim13 ∥ trim23 (nucmer) → MERGE_TRIMMED
 *                            → single CHECKV run → PICK_BEST_TRIM
 *   4. BLAST trim         (BLAST_TRIM subworkflow; controlled by params)
 *                            singletons + cluster-trim failures → BLAST → nucmer → CheckV
 *                            complete → recluster | no hit or incomplete → unvalidated
 *   5. Recluster          (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER → GET_CENTROIDS)
 *
 * Key params controlling step 4:
 *   --run_initial_blast   true/false  (default: true)
 *   --run_secondary_blast true/false  (default: false)
 *   --initial_blastdb     /path/to/db (required if run_initial_blast=true)
 *   --secondary_blastdb   /path/to/db (required if run_secondary_blast=true)
 *
 * Usage:
 *   nextflow run main.nf -profile slurm,conda [--fasta_list /path/to/fasta_list.txt] [--outdir /path]
 *   nextflow run main.nf -stub -profile local   # validate DAG
 */

nextflow.enable.dsl = 2

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
include { PIPELINE_SUMMARY                                  } from './modules/summary'

workflow {

    // -----------------------------------------------------------------------
    // Parameter validation
    // -----------------------------------------------------------------------
    def errors = []
    if (!params.fasta_list) errors << "  --fasta_list is required"
    if (!params.outdir)   errors << "  --outdir is required"
    if (!params.checkvdb) errors << "  --checkvdb is required"
    if (params.run_initial_blast && !params.initial_blastdb)
        errors << "  --initial_blastdb is required when --run_initial_blast=true"
    if (params.run_secondary_blast && !params.secondary_blastdb)
        errors << "  --secondary_blastdb is required when --run_secondary_blast=true"
    if (errors) {
        log.error "Parameter errors:\n" + errors.join("\n")
        System.exit(1)
    }

    // Validate BLAST db paths exist (skip during stub runs)
    if (!workflow.stubRun) {
        if (params.run_initial_blast && !file(params.initial_blastdb).exists())
            error "initial_blastdb not found: ${params.initial_blastdb}"
        if (params.run_secondary_blast && !file(params.secondary_blastdb).exists())
            error "secondary_blastdb not found: ${params.secondary_blastdb}"
    }

    log.info """
    ============================================
     wvdb_build: viral genome database pipeline
    ============================================
     fasta_list         : ${params.fasta_list}
     outdir             : ${params.outdir}
     ani                : ${params.ani}
     qcov               : ${params.qcov}
     completeness       : ${params.completeness}%
     threads            : ${params.threads}
     run_initial_blast  : ${params.run_initial_blast}
     run_secondary_blast: ${params.run_secondary_blast}
    ============================================
    """.stripIndent()

    def check = !workflow.stubRun
    checkvdb         = file(params.checkvdb, checkIfExists: check)
    fasta_list       = file(params.fasta_list, checkIfExists: check)
    // BLAST dbs passed as strings (not file objects) to prevent staging —
    // all index files (.nhr .nin .nsq etc.) must remain in the same directory
    initial_blastdb   = params.initial_blastdb   ?: ''
    secondary_blastdb = params.secondary_blastdb ?: ''

    // -----------------------------------------------------------------------
    // Step 1 — Collect and merge input genomes
    // -----------------------------------------------------------------------
    COLLECT_GENOMES(fasta_list)

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
    //   trim12, trim13, trim23 run in parallel → single CheckV run
    //   PICK_BEST_TRIM selects best per rank1 by completeness then length
    // -----------------------------------------------------------------------
    CLUSTER_TRIM(
        VCLUST_CLUSTER.out.clusters,
        COLLECT_GENOMES.out.lengths,
        COLLECT_GENOMES.out.merged_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Step 4 — BLAST-mode trimming
    //   Receives singletons + incomplete from branches A and B of step 3.
    //   Optional: controlled by run_initial_blast and run_secondary_blast.
    //   Sequences with no BLAST hit or incomplete after trimming → unvalidated.
    // -----------------------------------------------------------------------
    BLAST_TRIM(
        CLUSTER_TRIM.out.singletons,
        CLUSTER_TRIM.out.blast_trim_fasta,
        COLLECT_GENOMES.out.merged_fasta,
        initial_blastdb,
        secondary_blastdb,
        checkvdb,
        params.completeness
    )

    // -----------------------------------------------------------------------
    // Step 5 — Recluster complete representatives from steps 3 and 4
    //   Only trimmed + CheckV-complete sequences enter reclustering.
    //   Unvalidated sequences are written to unvalidated/ directory only.
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

    // -----------------------------------------------------------------------
    // Summary report — sequence counts at each step
    // -----------------------------------------------------------------------
    PIPELINE_SUMMARY(
        COLLECT_GENOMES.out.merged_fasta,
        VCLUST_CLUSTER.out.clusters,
        CLUSTER_TRIM.out.trimming_candidates,
        CLUSTER_TRIM.out.complete_fasta,
        CLUSTER_TRIM.out.trim_log,
        BLAST_TRIM.out.blast_trim_input,
        BLAST_TRIM.out.complete_fasta,
        BLAST_TRIM.out.unvalidated_fasta,
        BLAST_TRIM.out.unvalidated_report,
        recluster_input,
        GET_CENTROIDS.out.centroid_fasta
    )
}
