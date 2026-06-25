/*
 * subworkflows/blast_trim.nf
 * Step 4: BLAST-mode reference-guided trimming.
 *
 * Input: singletons + sequences incomplete after cluster_trim step.
 *
 * DAG:
 *                         ┌─ BLASTN (refseq_ev) ─ BLASTANI ─ TRIM_GENOMES_BLAST ─┐
 *   COLLECT_BLAST_INPUT ──┤                                                        ├─ MERGE_BLAST_TRIM ─ CHECKV ─ COMPLETENESS_FILTER
 *                         └─ BLASTN (metavr)    ─ BLASTANI ─ TRIM_GENOMES_BLAST ─┘
 *
 * The two BLAST database searches run in parallel.
 */

include { COLLECT_BLAST_INPUT                                  } from '../modules/blast_trim_tools'
include { BLASTN                                               } from '../modules/blast_trim_tools'
include { BLASTN              as BLASTN_METAVR                 } from '../modules/blast_trim_tools'
include { BLASTANI                                             } from '../modules/blast_trim_tools'
include { BLASTANI            as BLASTANI_METAVR               } from '../modules/blast_trim_tools'
include { TRIM_GENOMES_BLAST                                   } from '../modules/blast_trim_tools'
include { TRIM_GENOMES_BLAST  as TRIM_GENOMES_BLAST_METAVR     } from '../modules/blast_trim_tools'
include { MERGE_BLAST_TRIM                                     } from '../modules/blast_trim_tools'
include { CHECKV              } from '../modules/cluster_trim_tools'
include { COMPLETENESS_FILTER } from '../modules/cluster_trim_tools'

workflow BLAST_TRIM {

    take:
    singletons_ids        // path: singletons.txt from PARSE_CLUSTERS_ALL
    incomplete_fasta      // path: merged incomplete seqs from cluster_trim step
    all_fasta             // path: filtered_all.fasta
    refseq_ev_blastdb     // path: RefSeq + EsViritu BLAST db
    imgvr_blastdb         // path: IMG-VR v5 BLAST db
    checkvdb              // path: CheckV database
    completeness          // val:  completeness threshold (%)

    main:

    // -----------------------------------------------------------------------
    // Collect all blast_trim input sequences: singletons + cluster-trim failures
    // -----------------------------------------------------------------------
    COLLECT_BLAST_INPUT(
        singletons_ids,
        incomplete_fasta,
        all_fasta
    )

    // -----------------------------------------------------------------------
    // Parallel BLAST searches against refseq_ev and metavr databases
    // -----------------------------------------------------------------------
    BLASTN(
        COLLECT_BLAST_INPUT.out.blast_trim_fasta,
        refseq_ev_blastdb,
        "refseq_ev"
    )

    BLASTN_METAVR(
        COLLECT_BLAST_INPUT.out.blast_trim_fasta,
        imgvr_blastdb,
        "metavr"
    )

    // -----------------------------------------------------------------------
    // Compute ANI from BLAST output (parallel)
    // -----------------------------------------------------------------------
    BLASTANI(
        BLASTN.out.blastn_tsv
    )

    BLASTANI_METAVR(
        BLASTN_METAVR.out.blastn_tsv
    )

    // -----------------------------------------------------------------------
    // BLAST-mode trimming against each database (parallel)
    // -----------------------------------------------------------------------
    TRIM_GENOMES_BLAST(
        BLASTANI.out.ani_tsv,
        COLLECT_BLAST_INPUT.out.blast_trim_fasta,
        refseq_ev_blastdb,
        "refseq_ev"
    )

    TRIM_GENOMES_BLAST_METAVR(
        BLASTANI_METAVR.out.ani_tsv,
        COLLECT_BLAST_INPUT.out.blast_trim_fasta,
        imgvr_blastdb,
        "metavr"
    )

    // -----------------------------------------------------------------------
    // Merge: all refseq_ev trimmed + metavr-only additions
    // -----------------------------------------------------------------------
    MERGE_BLAST_TRIM(
        TRIM_GENOMES_BLAST.out.trimmed_fasta,
        TRIM_GENOMES_BLAST.out.trimming_bed,
        TRIM_GENOMES_BLAST_METAVR.out.trimmed_fasta,
        TRIM_GENOMES_BLAST_METAVR.out.trimming_bed
    )

    // -----------------------------------------------------------------------
    // CheckV quality assessment on merged trimmed set
    // -----------------------------------------------------------------------
    CHECKV(
        MERGE_BLAST_TRIM.out.trimmed_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Filter by completeness — incomplete seqs are not fed further
    // (no downstream step to fall back to after blast_trim)
    // -----------------------------------------------------------------------
    COMPLETENESS_FILTER(
        CHECKV.out.quality_summary,
        MERGE_BLAST_TRIM.out.trimmed_fasta,
        all_fasta,
        completeness,
        false
    )

    emit:
    complete_fasta = COMPLETENESS_FILTER.out.complete_fasta
    incomplete_ids = COMPLETENESS_FILTER.out.incomplete_ids
}
