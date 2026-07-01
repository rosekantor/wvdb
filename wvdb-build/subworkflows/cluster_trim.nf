/*
 * subworkflows/cluster_trim.nf
 * Step 3: Cluster-based trimming pipeline.
 *
 * DAG:
 *   PARSE_CLUSTERS_ALL
 *       ↓
 *   EXTRACT_TRIMMING_SEQS
 *       ↓
 *   MINIMAP2_PAIRWISE          ← all-vs-all on rank1/2/3 candidates
 *       ↓
 *   TRIM_FROM_PAF              ← picks best alignment per cluster across
 *                                 pairs12, pairs13, pairs23 by alignment
 *                                 block length; one trimmed.fasta output
 *       ↓
 *   CHECKV                     ← one run, no duplicate IDs
 *       ↓
 *   PICK_BEST_TRIM             ← complete → recluster | incomplete → blast_trim
 *       ↓
 *   FETCH_BLAST_INPUT_SEQS     ← retrieve untrimmed seqs for blast_trim step
 *
 * Emits:
 *   complete_fasta   — complete trimmed representatives  → step 5 (recluster)
 *   blast_trim_fasta — untrimmed incomplete sequences    → step 4 (blast_trim)
 *   singletons       — 1-member cluster IDs              → step 4 (blast_trim)
 *   paf_gz           — bgzipped PAF for network analysis
 */

include { PARSE_CLUSTERS_ALL    } from '../modules/cluster_trim_tools'
include { EXTRACT_TRIMMING_SEQS } from '../modules/cluster_trim_tools'
include { MINIMAP2_PAIRWISE     } from '../modules/cluster_trim_tools'
include { TRIM_FROM_PAF         } from '../modules/cluster_trim_tools'
include { CHECKV                } from '../modules/cluster_trim_tools'
include { PICK_BEST_TRIM        } from '../modules/cluster_trim_tools'
include { FETCH_BLAST_INPUT_SEQS } from '../modules/cluster_trim_tools'

workflow CLUSTER_TRIM {

    take:
    clusters      // path: vclust_clusters.tsv
    lengths       // path: filtered_all.length.txt
    all_fasta     // path: filtered_all.fasta
    checkvdb      // path: CheckV database

    main:

    // Parse clusters → all pair files + candidate IDs in one pass
    PARSE_CLUSTERS_ALL(
        clusters,
        lengths
    )

    // Extract rank1/2/3 sequences for minimap2
    EXTRACT_TRIMMING_SEQS(
        PARSE_CLUSTERS_ALL.out.candidate_ids,
        all_fasta
    )

    // All-vs-all minimap2 on trimming candidates
    MINIMAP2_PAIRWISE(
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta
    )

    // Pick best alignment per cluster across pairs12/13/23, write trimmed FASTA
    TRIM_FROM_PAF(
        MINIMAP2_PAIRWISE.out.paf,
        PARSE_CLUSTERS_ALL.out.pairs12,
        PARSE_CLUSTERS_ALL.out.pairs13,
        PARSE_CLUSTERS_ALL.out.pairs23,
        EXTRACT_TRIMMING_SEQS.out.candidates_fasta
    )

    // CheckV quality assessment on trimmed sequences
    CHECKV(
        TRIM_FROM_PAF.out.trimmed_fasta,
        checkvdb
    )

    // Split by completeness
    PICK_BEST_TRIM(
        TRIM_FROM_PAF.out.trimmed_fasta,
        CHECKV.out.quality_summary
    )

    // Retrieve untrimmed sequences for blast_trim step
    FETCH_BLAST_INPUT_SEQS(
        PICK_BEST_TRIM.out.blast_trim_ids,
        all_fasta
    )

    emit:
    complete_fasta   = PICK_BEST_TRIM.out.complete_fasta
    blast_trim_fasta = FETCH_BLAST_INPUT_SEQS.out.blast_trim_fasta
    singletons       = PARSE_CLUSTERS_ALL.out.singletons
    paf_gz           = MINIMAP2_PAIRWISE.out.paf_gz
}
