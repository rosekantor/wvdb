/*
 * subworkflows/characterize_proteins.nf
 * Protein characterization subworkflow.
 *
 * Input: predicted proteins .faa from geNomad
 * Steps:
 *   DIAMOND × N databases (parallel, if run_diamond=true)
 *   HMMSEARCH × N profiles  (parallel, if run_hmmsearch=true)
 *   PROTEIN_SUMMARY — per-contig hit count summary
 *
 * Databases and profiles are specified via CSV files:
 *   --diamond_dbs    diamond_dbs.csv   (name,path per row)
 *   --hmm_profiles   hmm_profiles.csv  (name,path per row)
 *
 * DAG:
 *   proteins_faa
 *       ├──► DIAMOND (initial_db) ──┐
 *       ├──► DIAMOND (nr)         ──┤
 *       ├──► ...                  ──┼──► PROTEIN_SUMMARY
 *       ├──► HMMSEARCH (pfam)     ──┤
 *       ├──► HMMSEARCH (vogdb)    ──┤
 *       └──► ...                  ──┘
 *
 * Emits:
 *   diamond_results_dir   — directory of *.diamond.tsv files
 *   hmmsearch_results_dir — directory of *.hmmsearch.tsv files
 *   protein_summary_tsv   — per-contig protein hit counts
 */

include { DIAMOND         } from '../modules/protein_tools'
include { HMMSEARCH       } from '../modules/protein_tools'
include { PROTEIN_SUMMARY } from '../modules/protein_tools'

workflow CHARACTERIZE_PROTEINS {

    take:
    proteins_faa          // path: geNomad predicted proteins .faa

    main:

    // -----------------------------------------------------------------------
    // DIAMOND searches (parallel across all databases in diamond_dbs.csv)
    // -----------------------------------------------------------------------
    if ( params.run_diamond ) {
        // Filter unified databases CSV to diamond rows
        diamond_dbs_ch = Channel
            .fromPath(params.databases)
            .splitCsv(header: true)
            .filter { row -> row.type == "diamond" }
            .map { row -> tuple(row.name, row.path) }

        DIAMOND(
            proteins_faa,
            diamond_dbs_ch.map { name, path -> name },
            diamond_dbs_ch.map { name, path -> path }
        )

        diamond_tsvs = DIAMOND.out.diamond_result
            .map { db_name, tsv -> tsv }
            .collect()
            .ifEmpty([file('NO_FILE_DIAMOND')])

    } else {
        diamond_tsvs = Channel.value([file('NO_FILE_DIAMOND')])
    }

    // -----------------------------------------------------------------------
    // HMMSEARCH against all profiles in hmm_profiles.csv (parallel)
    // -----------------------------------------------------------------------
    if ( params.run_hmmsearch ) {
        // Filter unified databases CSV to hmm rows
        hmm_profiles_ch = Channel
            .fromPath(params.databases)
            .splitCsv(header: true)
            .filter { row -> row.type == "hmm" }
            .map { row -> tuple(row.name, file(row.path)) }

        HMMSEARCH(
            proteins_faa,
            hmm_profiles_ch.map { name, path -> name },
            hmm_profiles_ch.map { name, path -> path }
        )

        hmmsearch_tsvs = HMMSEARCH.out.hmmsearch_result
            .map { profile_name, tsv -> tsv }
            .collect()
            .ifEmpty([file('NO_FILE_HMMSEARCH')])

    } else {
        hmmsearch_tsvs = Channel.value([file('NO_FILE_HMMSEARCH')])
    }

    // -----------------------------------------------------------------------
    // Per-contig protein hit summary
    // -----------------------------------------------------------------------
    PROTEIN_SUMMARY(
        proteins_faa,
        diamond_tsvs,
        hmmsearch_tsvs
    )

    emit:
    diamond_tsvs         = diamond_tsvs
    hmmsearch_tsvs       = hmmsearch_tsvs
    protein_summary_tsv  = PROTEIN_SUMMARY.out.protein_summary_tsv
}
