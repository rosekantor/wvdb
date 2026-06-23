/*
 * subworkflows/branchC.nf
 * Branch C: BLAST-mode trimming of singletons and cluster reps that were
 * incomplete after Branches A and B.
 *
 * DAG:
 *                         ┌─ BLASTN (refseq_ev) ─ BLASTANI ─ TRIM_GENOMES_BLAST ─┐
 *   COLLECT_BRANCHC ──────┤                                                        ├─ MERGE_BRANCHC ─ CHECKV ─ COMPLETENESS_FILTER
 *                         └─ BLASTN (metavr)    ─ BLASTANI ─ TRIM_GENOMES_BLAST ─┘
 *
 * The two BLAST branches are independent and run in parallel.
 */

include { COLLECT_BRANCHC    } from '../modules/blast_trim'
include { BLASTN              } from '../modules/blast_trim'
include { BLASTANI            } from '../modules/blast_trim'
include { TRIM_GENOMES_BLAST  } from '../modules/blast_trim'
include { BLASTN              as BLASTN_METAVR             } from '../modules/blast_trim'
include { BLASTANI            as BLASTANI_METAVR           } from '../modules/blast_trim'
include { TRIM_GENOMES_BLAST  as TRIM_GENOMES_BLAST_METAVR } from '../modules/blast_trim'
include { MERGE_BRANCHC       } from '../modules/blast_trim'
include { CHECKV              } from '../modules/trim_filter'
include { COMPLETENESS_FILTER } from '../modules/trim_filter'

workflow BRANCH_C {

    take:
    singletons_ids        // path: cluster_singletons.txt from BRANCH_A
    incomplete_fasta      // path: cluster_reps_notcomplete_after23.fasta from BRANCH_B
    all_fasta             // path: filtered_all.fasta
    refseq_ev_blastdb     // path: RefSeq + EsViritu BLAST db
    imgvr_blastdb         // path: IMG-VR v5 BLAST db
    checkvdb              // path: CheckV database
    completeness          // val:  completeness threshold (%)

    main:

    // -----------------------------------------------------------------------
    // Collect branch C input genomes
    // -----------------------------------------------------------------------
    COLLECT_BRANCHC(
        singletons_ids,
        incomplete_fasta,
        all_fasta
    )

    // -----------------------------------------------------------------------
    // Parallel BLAST searches against refseq_ev and metavr databases
    // -----------------------------------------------------------------------
    BLASTN(
        COLLECT_BRANCHC.out.branchC_fasta,
        refseq_ev_blastdb,
        "refseq_ev"
    )

    BLASTN_METAVR(
        COLLECT_BRANCHC.out.branchC_fasta,
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
    // Blast-mode trimming against each database (parallel)
    // -----------------------------------------------------------------------
    TRIM_GENOMES_BLAST(
        BLASTANI.out.ani_tsv,
        COLLECT_BRANCHC.out.branchC_fasta,
        refseq_ev_blastdb,
        "refseq_ev"
    )

    TRIM_GENOMES_BLAST_METAVR(
        BLASTANI_METAVR.out.ani_tsv,
        COLLECT_BRANCHC.out.branchC_fasta,
        imgvr_blastdb,
        "metavr"
    )

    // -----------------------------------------------------------------------
    // Merge: keep refseq_ev trimmed + metavr-only additions
    // -----------------------------------------------------------------------
    MERGE_BRANCHC(
        TRIM_GENOMES_BLAST.out.trimmed_fasta,
        TRIM_GENOMES_BLAST.out.trimming_bed,
        TRIM_GENOMES_BLAST_METAVR.out.trimmed_fasta,
        TRIM_GENOMES_BLAST_METAVR.out.trimming_bed
    )

    // -----------------------------------------------------------------------
    // CheckV quality assessment on merged trimmed set
    // -----------------------------------------------------------------------
    CHECKV(
        MERGE_BRANCHC.out.trimmed_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Filter by completeness
    // -----------------------------------------------------------------------
    COMPLETENESS_FILTER(
        CHECKV.out.quality_summary,
        MERGE_BRANCHC.out.trimmed_fasta,
        all_fasta,
        completeness,
        false              // branch C does not feed another branch; no incomplete FASTA needed
    )

    emit:
    complete_fasta = COMPLETENESS_FILTER.out.complete_fasta
    incomplete_ids = COMPLETENESS_FILTER.out.incomplete_ids
}
