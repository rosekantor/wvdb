/*
 * subworkflows/cluster_trim.nf
 * Step 3: Cluster-based trimming pipeline.
 *
 * Three sequential branches:
 *
 *   Branch A (parallel trim12 + trim13):
 *     PARSE_CLUSTERS_ALL → EXTRACT_TRIMMING_SEQS
 *         ↓──────────────────┐
 *     TRIM_GENOMES (trim12)  TRIM_GENOMES_13 (trim13)   ← parallel
 *         └──────────────────┘
 *         PREPARE_CHECKV_INPUT → CHECKV → PICK_BEST_TRIM1
 *
 *   Branch B (sequential, incomplete rank1 subset only):
 *     FILTER_PAIRS23 → TRIM_GENOMES_23 → CHECKV_TRIM23 → COMPLETENESS_FILTER_TRIM23
 *
 *   Outputs to step 4 (blast_trim):
 *     FETCH_BLAST_INPUT_SEQS  — 2-member cluster trim1 failures
 *     COMPLETENESS_FILTER_TRIM23.incomplete_fasta  — trim23 failures
 *     singletons              — 1-member clusters
 *
 * Emits:
 *   complete_fasta      — complete reps from branch A     → step 5 (recluster)
 *   complete_fasta23    — complete reps from branch B     → step 5 (recluster)
 *   blast_trim_fasta    — untrimmed seqs from 2-member failures → step 4
 *   blast_trim_fasta23  — untrimmed seqs from branch B failures → step 4
 *   singletons          — 1-member cluster IDs             → step 4
 */

include { PARSE_CLUSTERS_ALL                                   } from '../modules/cluster_trim_tools'
include { EXTRACT_TRIMMING_SEQS                                } from '../modules/cluster_trim_tools'
include { TRIM_GENOMES                                         } from '../modules/cluster_trim_tools'
include { TRIM_GENOMES    as TRIM_GENOMES_13                   } from '../modules/cluster_trim_tools'
include { TRIM_GENOMES    as TRIM_GENOMES_23                   } from '../modules/cluster_trim_tools'
include { PREPARE_CHECKV_INPUT                                 } from '../modules/cluster_trim_tools'
include { PICK_BEST_TRIM1                                      } from '../modules/cluster_trim_tools'
include { FETCH_BLAST_INPUT_SEQS                               } from '../modules/cluster_trim_tools'
include { FILTER_PAIRS23                                       } from '../modules/cluster_trim_tools'
include { CHECKV                                               } from '../modules/cluster_trim_tools'
include { CHECKV          as CHECKV_TRIM23                     } from '../modules/cluster_trim_tools'
include { COMPLETENESS_FILTER                                  } from '../modules/cluster_trim_tools'
include { COMPLETENESS_FILTER as COMPLETENESS_FILTER_TRIM23    } from '../modules/cluster_trim_tools'

workflow CLUSTER_TRIM {

    take:
    clusters      // path: vclust_clusters.tsv
    lengths       // path: filtered_all.length.txt
    all_fasta     // path: filtered_all.fasta
    checkvdb      // path: CheckV database

    main:

    // -----------------------------------------------------------------------
    // Shared setup: parse clusters + extract trimming candidates
    // -----------------------------------------------------------------------
    PARSE_CLUSTERS_ALL(
        clusters,
        lengths
    )

    EXTRACT_TRIMMING_SEQS(
        PARSE_CLUSTERS_ALL.out.candidate_ids,
        all_fasta
    )

    // -----------------------------------------------------------------------
    // Branch A: trim12 and trim13 run in parallel
    // Both use the pre-indexed trimming_candidates.fasta via trim_genomes.py
    // -----------------------------------------------------------------------
    TRIM_GENOMES(
        PARSE_CLUSTERS_ALL.out.pairs12,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta,
        "trim12"
    )

    TRIM_GENOMES_13(
        PARSE_CLUSTERS_ALL.out.pairs13,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta,
        "trim13"
    )

    // Prepare combined FASTA with suffixed IDs for CheckV
    // Both trim12 and trim13 assessed together so completeness drives selection
    PREPARE_CHECKV_INPUT(
        TRIM_GENOMES.out.trimmed_fasta,
        TRIM_GENOMES_13.out.trimmed_fasta
    )

    // CheckV on combined trim1 candidates (<rank1_id>_trim12 and _trim13)
    CHECKV(
        PREPARE_CHECKV_INPUT.out.checkv_input_fasta,
        checkvdb
    )

    // Split by completeness:
    //   complete            → step 5 (recluster)
    //   incomplete_ids      → branch B (trim23)
    //   blast_trim_ids      → step 4 (blast_trim); 2-member clusters with no rank3
    PICK_BEST_TRIM1(
        TRIM_GENOMES.out.trimmed_fasta,
        TRIM_GENOMES_13.out.trimmed_fasta,
        CHECKV.out.quality_summary,
        PARSE_CLUSTERS_ALL.out.pairs12,
        PARSE_CLUSTERS_ALL.out.pairs13
    )

    // Retrieve untrimmed seqs for 2-member cluster failures → step 4
    FETCH_BLAST_INPUT_SEQS(
        PICK_BEST_TRIM1.out.blast_trim_ids,
        all_fasta
    )

    // -----------------------------------------------------------------------
    // Branch B: trim23 runs sequentially, only on incomplete rank1 clusters
    // FILTER_PAIRS23 restricts pairs23.tsv before nucmer runs, saving compute
    // -----------------------------------------------------------------------
    FILTER_PAIRS23(
        PARSE_CLUSTERS_ALL.out.pairs23,
        PARSE_CLUSTERS_ALL.out.pairs12,
        PICK_BEST_TRIM1.out.incomplete_ids
    )

    TRIM_GENOMES_23(
        FILTER_PAIRS23.out.pairs23_filtered,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta,
        "trim23"
    )

    CHECKV_TRIM23(
        TRIM_GENOMES_23.out.trimmed_fasta,
        checkvdb
    )

    // Split trim23 by completeness:
    //   complete            → step 5 (recluster)
    //   incomplete_fasta    → step 4 (blast_trim)
    COMPLETENESS_FILTER_TRIM23(
        CHECKV_TRIM23.out.quality_summary,
        TRIM_GENOMES_23.out.trimmed_fasta,
        all_fasta,
        params.completeness,
        true    // emit untrimmed incomplete seqs for blast_trim step
    )

    emit:
    complete_fasta      = PICK_BEST_TRIM1.out.complete_fasta
    complete_fasta23    = COMPLETENESS_FILTER_TRIM23.out.complete_fasta
    blast_trim_fasta    = FETCH_BLAST_INPUT_SEQS.out.blast_trim_fasta
    blast_trim_fasta23  = COMPLETENESS_FILTER_TRIM23.out.incomplete_fasta
    singletons          = PARSE_CLUSTERS_ALL.out.singletons
}
