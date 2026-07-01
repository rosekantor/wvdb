/*
 * subworkflows/cluster_trim.nf
 * Cluster-based trimming pipeline (steps 3/4), replacing the old
 * BRANCH_A / BRANCH_B subworkflows.
 *
 * DAG:
 *   PARSE_CLUSTERS_ALL
 *       ↓
 *   EXTRACT_TRIMMING_SEQS
 *       ↓
 *   MINIMAP2_PAIRWISE
 *       ↓──────────────┬──────────────┐
 *   TRIM_GENOMES_12  TRIM_GENOMES_13  TRIM_GENOMES_23   (parallel)
 *       ↓              ↓
 *   MERGE_BEST_TRIM1   (longer of trim12 vs trim13 per rank1)
 *       ↓
 *   CHECKV_TRIM1
 *       ↓
 *   PICK_BEST_TRIM1    (split: complete / incomplete-try-23 / incomplete-branch-C)
 *       ↓
 *   FILTER_TRIM23      (restrict trim23 to incomplete rank1 clusters)
 *       ↓
 *   CHECKV_TRIM23
 *       ↓
 *   COMPLETENESS_FILTER_TRIM23
 *
 * Emits:
 *   complete_fasta      — complete reps from trim1 step       → step 6
 *   complete_fasta23    — complete reps from trim23 step      → step 6
 *   branch_c_fasta      — untrimmed seqs for branch C (trim1 failures)
 *   branch_c_fasta23    — untrimmed seqs for branch C (trim23 failures)
 *   singletons          — singleton IDs                       → branch C
 *   paf_gz              — full bgzipped PAF for network analysis
 */

include { PARSE_CLUSTERS_ALL                          } from '../modules/trim_filter'
include { EXTRACT_TRIMMING_SEQS                       } from '../modules/trim_filter'
include { MINIMAP2_PAIRWISE                           } from '../modules/trim_filter'
include { TRIM_GENOMES                                } from '../modules/trim_filter'
include { TRIM_GENOMES as TRIM_GENOMES_13             } from '../modules/trim_filter'
include { TRIM_GENOMES as TRIM_GENOMES_23             } from '../modules/trim_filter'
include { MERGE_BEST_TRIM1                            } from '../modules/trim_filter'
include { PICK_BEST_TRIM1                             } from '../modules/trim_filter'
include { FILTER_TRIM23                               } from '../modules/trim_filter'
include { FETCH_BRANCH_C_SEQS                         } from '../modules/trim_filter'
include { CHECKV                                      } from '../modules/trim_filter'
include { CHECKV          as CHECKV_TRIM23            } from '../modules/trim_filter'
include { COMPLETENESS_FILTER                         } from '../modules/trim_filter'
include { COMPLETENESS_FILTER as COMPLETENESS_FILTER_TRIM23 } from '../modules/trim_filter'

workflow CLUSTER_TRIM {

    take:
    clusters      // path: vclust_clusters.tsv
    lengths       // path: filtered_all.length.txt
    all_fasta     // path: filtered_all.fasta
    checkvdb      // path: CheckV database

    main:

    // -----------------------------------------------------------------------
    // Parse clusters → all pair files + candidate IDs in one pass
    // -----------------------------------------------------------------------
    PARSE_CLUSTERS_ALL(
        clusters,
        lengths
    )

    // -----------------------------------------------------------------------
    // Extract rank1/2/3 sequences for minimap2
    // -----------------------------------------------------------------------
    EXTRACT_TRIMMING_SEQS(
        PARSE_CLUSTERS_ALL.out.candidate_ids,
        all_fasta
    )

    // -----------------------------------------------------------------------
    // Single minimap2 all-vs-all on trimming candidates
    // -----------------------------------------------------------------------
    MINIMAP2_PAIRWISE(
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta
    )

    // -----------------------------------------------------------------------
    // Parallel trimming: trim12, trim13, trim23 all filter the same PAF
    // -----------------------------------------------------------------------
    TRIM_GENOMES(
        MINIMAP2_PAIRWISE.out.paf,
        PARSE_CLUSTERS_ALL.out.pairs12,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta,
        "trim12"
    )

    TRIM_GENOMES_13(
        MINIMAP2_PAIRWISE.out.paf,
        PARSE_CLUSTERS_ALL.out.pairs13,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta,
        "trim13"
    )

    TRIM_GENOMES_23(
        MINIMAP2_PAIRWISE.out.paf,
        PARSE_CLUSTERS_ALL.out.pairs23,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta,
        "trim23"
    )

    // -----------------------------------------------------------------------
    // Merge trim12 + trim13; keep longer per rank1 ID
    // -----------------------------------------------------------------------
    MERGE_BEST_TRIM1(
        TRIM_GENOMES.out.trimmed_fasta,
        TRIM_GENOMES_13.out.trimmed_fasta
    )

    // -----------------------------------------------------------------------
    // CheckV on merged trim1 candidates
    // -----------------------------------------------------------------------
    CHECKV(
        MERGE_BEST_TRIM1.out.best_trim1_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Pick best trim1 per cluster; split by completeness
    // -----------------------------------------------------------------------
    PICK_BEST_TRIM1(
        TRIM_GENOMES.out.trimmed_fasta,
        TRIM_GENOMES_13.out.trimmed_fasta,
        CHECKV.out.quality_summary,
        PARSE_CLUSTERS_ALL.out.pairs12,
        PARSE_CLUSTERS_ALL.out.pairs13
    )

    // -----------------------------------------------------------------------
    // Retrieve untrimmed seqs for 2-member clusters that failed trim1 → branch C
    // -----------------------------------------------------------------------
    FETCH_BRANCH_C_SEQS(
        PICK_BEST_TRIM1.out.branch_c_ids,
        all_fasta
    )

    // -----------------------------------------------------------------------
    // Filter trim23 to only clusters where trim1 was incomplete
    // -----------------------------------------------------------------------
    FILTER_TRIM23(
        TRIM_GENOMES_23.out.trimmed_fasta,
        PICK_BEST_TRIM1.out.incomplete_ids,
        PARSE_CLUSTERS_ALL.out.pairs23
    )

    // -----------------------------------------------------------------------
    // CheckV on filtered trim23
    // -----------------------------------------------------------------------
    CHECKV_TRIM23(
        FILTER_TRIM23.out.filtered_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Completeness filter on trim23 results
    // emit_incomplete_fasta=true: writes untrimmed seqs for branch C
    // -----------------------------------------------------------------------
    COMPLETENESS_FILTER_TRIM23(
        CHECKV_TRIM23.out.quality_summary,
        FILTER_TRIM23.out.filtered_fasta,
        all_fasta,
        params.completeness,
        true
    )

    emit:
    complete_fasta    = PICK_BEST_TRIM1.out.complete_fasta
    complete_fasta23  = COMPLETENESS_FILTER_TRIM23.out.complete_fasta
    branch_c_fasta    = FETCH_BRANCH_C_SEQS.out.branch_c_fasta
    branch_c_fasta23  = COMPLETENESS_FILTER_TRIM23.out.incomplete_fasta
    singletons        = PARSE_CLUSTERS_ALL.out.singletons
    paf_gz            = MINIMAP2_PAIRWISE.out.paf_gz
}
