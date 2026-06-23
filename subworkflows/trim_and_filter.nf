/*
 * subworkflows/trim_and_filter.nf
 * Shared subworkflow called for both Branch A (rank12) and Branch B (rank23).
 *
 *   PARSE_CLUSTERS → TRIM_GENOMES → CHECKV → COMPLETENESS_FILTER
 *
 * task.ext.publish_dir is set via nextflow.config withName selectors,
 * scoped to BRANCH_A:* and BRANCH_B:* so each process publishes to the
 * correct subdirectory (3_branchA or 4_branchB) without needing a val input.
 */

include { PARSE_CLUSTERS      } from '../modules/trim_filter'
include { TRIM_GENOMES        } from '../modules/trim_filter'
include { CHECKV              } from '../modules/trim_filter'
include { COMPLETENESS_FILTER } from '../modules/trim_filter'

workflow TRIM_AND_FILTER {

    take:
    clusters              // path: vclust_clusters.tsv
    lengths               // path: filtered_all.length.txt
    all_fasta             // path: filtered_all.fasta
    checkvdb              // path: checkv database dir
    mode                  // val:  "rank12" or "rank23"
    restrict_ids          // path: incomplete IDs to restrict to, or file('NO_FILE')
    completeness          // val:  completeness threshold (%)
    emit_incomplete_fasta // val:  true for branch B, false for branch A

    main:
    PARSE_CLUSTERS(
        clusters,
        lengths,
        mode,
        restrict_ids
    )
    PARSE_CLUSTERS.out.pairs.dump(tag: 'parse_clusters_pairs')

    TRIM_GENOMES(
        PARSE_CLUSTERS.out.pairs,
        all_fasta
    )

    CHECKV(
        TRIM_GENOMES.out.trimmed_fasta,
        checkvdb
    )

    COMPLETENESS_FILTER(
        CHECKV.out.quality_summary,
        TRIM_GENOMES.out.trimmed_fasta,
        all_fasta,
        completeness,
        emit_incomplete_fasta
    )

    emit:
    complete_fasta   = COMPLETENESS_FILTER.out.complete_fasta
    incomplete_ids   = COMPLETENESS_FILTER.out.incomplete_ids
    incomplete_fasta = COMPLETENESS_FILTER.out.incomplete_fasta
    singletons       = PARSE_CLUSTERS.out.singletons
}
