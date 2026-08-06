/*
 * subworkflows/blast_trim.nf
 * Step 4: BLAST-mode reference-guided trimming.
 *
 * Controlled by params:
 *   params.run_initial_blast   (default: true)  — search initial_blastdb
 *   params.run_secondary_blast (default: false) — search secondary_blastdb
 *
 * DAG (both blast steps enabled):
 *   COLLECT_BLAST_INPUT → GET_QUERY_IDS
 *       ├──► BLASTN (initial)   → BLASTANI_INITIAL   ──┐
 *       └──► BLASTN (secondary) → BLASTANI_SECONDARY ──┴──► SELECT_BEST_BLAST_HIT
 *                                                                    │
 *                                        ┌───────────────────────────┤
 *                                        ▼                           ▼                ▼
 *                                 initial_hits.tsv     secondary_hits.tsv      no_hit_ids.txt
 *                                        │                           │                │
 *                                TRIM_GENOMES_BLAST   TRIM_GENOMES_BLAST_SECONDARY    │
 *                                        │                           │                │
 *                                        └───────────┬───────────────┘                │
 *                                                    ▼                                │
 *                                             MERGE_TRIMMED                           │
 *                                                    ▼                                │
 *                                                  CHECKV                             │
 *                                                    ▼                                │
 *                                        COMPLETENESS_FILTER                          │
 *                                         ┌──────┴──────┐                            │
 *                                         ▼             ▼                            │
 *                                      complete    incomplete ──────────────────────► │
 *                                         │                                           │
 *                                         ▼                                           ▼
 *                                     recluster                            COLLECT_UNVALIDATED
 *                                                                           → unvalidated/
 *
 * When run_initial_blast=false:
 *   All inputs go directly to COLLECT_UNVALIDATED (no trimming possible).
 *
 * When run_secondary_blast=false (default):
 *   Secondary BLAST step skipped; secondary_hits.tsv is empty.
 */

include { COLLECT_BLAST_INPUT                                      } from '../modules/blast_trim_tools'
include { GET_QUERY_IDS                                            } from '../modules/blast_trim_tools'
include { BLASTN                                                   } from '../modules/blast_trim_tools'
include { BLASTN              as BLASTN_SECONDARY                  } from '../modules/blast_trim_tools'
include { BLASTANI                                                 } from '../modules/blast_trim_tools'
include { BLASTANI            as BLASTANI_SECONDARY                } from '../modules/blast_trim_tools'
include { SELECT_BEST_BLAST_HIT                                    } from '../modules/blast_trim_tools'
include { TRIM_GENOMES_BLAST                                       } from '../modules/blast_trim_tools'
include { TRIM_GENOMES_BLAST  as TRIM_GENOMES_BLAST_SECONDARY      } from '../modules/blast_trim_tools'
include { MERGE_TRIMMED                                            } from '../modules/blast_trim_tools'
include { COLLECT_UNVALIDATED                                      } from '../modules/blast_trim_tools'
include { CHECKV              } from '../modules/cluster_trim_tools'
include { COMPLETENESS_FILTER } from '../modules/cluster_trim_tools'

workflow BLAST_TRIM {

    take:
    singletons_ids        // path: singletons.txt from PARSE_CLUSTERS_ALL
    incomplete_fasta      // path: merged incomplete seqs from cluster_trim step
    all_fasta             // path: filtered_all.fasta
    initial_blastdb       // val:  initial BLAST db path (null if run_initial_blast=false)
    secondary_blastdb     // val:  secondary BLAST db path (null if run_secondary_blast=false)
    checkvdb              // path: CheckV database
    completeness          // val:  completeness threshold (%)

    main:

    // -----------------------------------------------------------------------
    // Collect blast_trim input sequences regardless of whether blast runs
    // -----------------------------------------------------------------------
    COLLECT_BLAST_INPUT(
        singletons_ids,
        incomplete_fasta,
        all_fasta
    )

    GET_QUERY_IDS(
        COLLECT_BLAST_INPUT.out.blast_trim_fasta
    )

    if ( params.run_initial_blast ) {

        // -------------------------------------------------------------------
        // BLASTN searches — both run in parallel if secondary enabled
        // -------------------------------------------------------------------
        BLASTN(
            COLLECT_BLAST_INPUT.out.blast_trim_fasta,
            initial_blastdb,
            "initial"
        )

        BLASTANI(
            BLASTN.out.blastn_tsv
        )

        // Secondary BLAST: only if run_secondary_blast=true
        if ( params.run_secondary_blast ) {
            BLASTN_SECONDARY(
                COLLECT_BLAST_INPUT.out.blast_trim_fasta,
                secondary_blastdb,
                "secondary"
            )
            BLASTANI_SECONDARY(
                BLASTN_SECONDARY.out.blastn_tsv
            )
            secondary_ani = BLASTANI_SECONDARY.out.ani_tsv
        } else {
            // Emit an empty file so SELECT_BEST_BLAST_HIT receives a valid path
            secondary_ani = Channel.of('no_secondary')
                .map { _x ->
                    def f = file("${workDir}/empty_secondary_ani.tsv")
                    f.text = ''
                    return f
                }
        }

        // -------------------------------------------------------------------
        // Route queries: initial preferred, secondary if no initial hit
        // -------------------------------------------------------------------
        SELECT_BEST_BLAST_HIT(
            BLASTANI.out.ani_tsv,
            secondary_ani,
            GET_QUERY_IDS.out.query_ids
        )

        // -------------------------------------------------------------------
        // Trimming: initial runs first; secondary sequential after
        // (secondary input depends on SELECT_BEST_BLAST_HIT output)
        // -------------------------------------------------------------------
        TRIM_GENOMES_BLAST(
            SELECT_BEST_BLAST_HIT.out.initial_hits,
            COLLECT_BLAST_INPUT.out.blast_trim_fasta,
            initial_blastdb,
            "initial"
        )

        if ( params.run_secondary_blast ) {
            TRIM_GENOMES_BLAST_SECONDARY(
                SELECT_BEST_BLAST_HIT.out.secondary_hits,
                COLLECT_BLAST_INPUT.out.blast_trim_fasta,
                secondary_blastdb,
                "secondary"
            )
            secondary_trimmed = TRIM_GENOMES_BLAST_SECONDARY.out.trimmed_fasta
            secondary_no_aln  = TRIM_GENOMES_BLAST_SECONDARY.out.no_alignment_ids
        } else {
            secondary_trimmed = Channel.of('no_secondary')
                .map { _x ->
                    def f = file("${workDir}/empty_secondary.fasta")
                    f.text = ''
                    return f
                }
            secondary_no_aln = Channel.of('no_secondary')
                .map { _x ->
                    def f = file("${workDir}/empty_secondary_no_aln.txt")
                    f.text = ''
                    return f
                }
        }

        MERGE_TRIMMED(
            TRIM_GENOMES_BLAST.out.trimmed_fasta,
            secondary_trimmed
        )

        CHECKV(
            MERGE_TRIMMED.out.trimmed_fasta,
            checkvdb
        )

        COMPLETENESS_FILTER(
            CHECKV.out.quality_summary,
            MERGE_TRIMMED.out.trimmed_fasta,
            all_fasta,
            completeness,
            false
        )

        complete_fasta   = COMPLETENESS_FILTER.out.complete_fasta
        incomplete_ids_ch = COMPLETENESS_FILTER.out.incomplete_ids
        no_hit_ids_ch    = SELECT_BEST_BLAST_HIT.out.no_hit_ids

        // Merge no-usable-nucmer-alignment IDs from initial + secondary
        no_alignment_ids_ch = TRIM_GENOMES_BLAST.out.no_alignment_ids
            .mix(secondary_no_aln)
            .collectFile(name: 'no_nucmer_alignment_merged.txt')

    } else {
        // run_initial_blast = false: no trimming, all sequences → unvalidated
        complete_fasta    = Channel.empty()
        incomplete_ids_ch = Channel.empty()
        no_hit_ids_ch     = GET_QUERY_IDS.out.query_ids
        no_alignment_ids_ch = Channel.of('no_blast').map { _x ->
            def f = file("${workDir}/empty_no_alignment.txt"); f.text = ''; return f }
    }

    // -----------------------------------------------------------------------
    // Collect all unvalidated sequences (no qualifying hit + no usable
    // nucmer alignment + incomplete after trim)
    // -----------------------------------------------------------------------
    def incomplete_ids_final = (params.run_initial_blast)
        ? incomplete_ids_ch
        : Channel.of('no_blast').map { _x ->
            def f = file("${workDir}/empty_incomplete.txt"); f.text = ''; return f }

    COLLECT_UNVALIDATED(
        no_hit_ids_ch,
        no_alignment_ids_ch,
        incomplete_ids_final,
        all_fasta
    )

    emit:
    complete_fasta       = complete_fasta
    blast_trim_input     = COLLECT_BLAST_INPUT.out.blast_trim_fasta
    unvalidated_fasta    = COLLECT_UNVALIDATED.out.unvalidated_fasta
    unvalidated_report   = COLLECT_UNVALIDATED.out.unvalidated_report
}
