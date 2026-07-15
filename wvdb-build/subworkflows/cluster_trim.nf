/*
 * subworkflows/cluster_trim.nf
 * Step 3: Cluster-based trimming pipeline.
 *
 * All three trim passes run in parallel from trimming_candidates.fasta.
 * A single CheckV run covers all trimmings via suffixed IDs.
 *
 * DAG:
 *   PARSE_CLUSTERS_ALL → EXTRACT_TRIMMING_SEQS
 *       ├──► TRIM_GENOMES    (trim12, rank1 vs rank2)  ──┐
 *       ├──► TRIM_GENOMES_13 (trim13, rank1 vs rank3)  ──┼──► MERGE_TRIMMED
 *       └──► TRIM_GENOMES_23 (trim23, rank2 vs rank3)  ──┘         │
 *                                                              CHECKV (one run)
 *                                                                    │
 *                                                            PICK_BEST_TRIM
 *                                                            (best of 3 per rank1)
 *                                                             ┌──────┴──────┐
 *                                                          complete      incomplete
 *                                                             │               │
 *                                                          step 5     FETCH_BLAST_INPUT_SEQS
 *                                                       (recluster)      → step 4
 *
 * Emits:
 *   complete_fasta       — complete reps → step 5 (recluster)
 *   blast_trim_fasta     — untrimmed incomplete seqs → step 4 (blast_trim)
 *   singletons           — 1-member cluster IDs → step 4 (blast_trim)
 *   trimming_candidates  — candidates.fasta (for PIPELINE_SUMMARY)
 */

include { PARSE_CLUSTERS_ALL    } from '../modules/cluster_trim_tools'
include { EXTRACT_TRIMMING_SEQS } from '../modules/cluster_trim_tools'
include { TRIM_GENOMES          } from '../modules/cluster_trim_tools'
include { TRIM_GENOMES          as TRIM_GENOMES_13 } from '../modules/cluster_trim_tools'
include { TRIM_GENOMES          as TRIM_GENOMES_23 } from '../modules/cluster_trim_tools'
include { MERGE_TRIMMED         } from '../modules/cluster_trim_tools'
include { CHECKV                } from '../modules/cluster_trim_tools'
include { PICK_BEST_TRIM        } from '../modules/cluster_trim_tools'
include { FETCH_BLAST_INPUT_SEQS} from '../modules/cluster_trim_tools'

workflow CLUSTER_TRIM {

    take:
    clusters      // path: vclust_clusters.tsv
    lengths       // path: filtered_all.length.txt
    all_fasta     // path: filtered_all.fasta
    checkvdb      // path: CheckV database

    main:

    // -----------------------------------------------------------------------
    // Parse clusters and extract trimming candidates
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
    // Three trim passes in parallel — all use pre-indexed candidates.fasta
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

    TRIM_GENOMES_23(
        PARSE_CLUSTERS_ALL.out.pairs23,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta,
        "trim23"
    )

    // -----------------------------------------------------------------------
    // Merge all three trim outputs with suffixed IDs → single CheckV run
    // -----------------------------------------------------------------------
    MERGE_TRIMMED(
        TRIM_GENOMES.out.trimmed_fasta,
        TRIM_GENOMES_13.out.trimmed_fasta,
        TRIM_GENOMES_23.out.trimmed_fasta,
        PARSE_CLUSTERS_ALL.out.pairs12
    )

    CHECKV(
        MERGE_TRIMMED.out.all_trimmed_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Select best trimming per rank1 by CheckV completeness then length
    // -----------------------------------------------------------------------
    PICK_BEST_TRIM(
        TRIM_GENOMES.out.trimmed_fasta,
        TRIM_GENOMES_13.out.trimmed_fasta,
        TRIM_GENOMES_23.out.trimmed_fasta,
        CHECKV.out.quality_summary,
        PARSE_CLUSTERS_ALL.out.pairs12
    )

    // Retrieve untrimmed sequences for incomplete rank1 clusters → step 4
    FETCH_BLAST_INPUT_SEQS(
        PICK_BEST_TRIM.out.blast_trim_ids,
        all_fasta
    )

    emit:
    complete_fasta      = PICK_BEST_TRIM.out.complete_fasta
    blast_trim_fasta    = FETCH_BLAST_INPUT_SEQS.out.blast_trim_fasta
    singletons          = PARSE_CLUSTERS_ALL.out.singletons
    trimming_candidates = EXTRACT_TRIMMING_SEQS.out.candidates_fasta
    trim_log            = PICK_BEST_TRIM.out.trim_log
}
