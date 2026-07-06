#!/usr/bin/env nextflow
/*
 * main.nf  —  wvdb_annotate: viral genome annotation pipeline
 *
 * Steps (all run in parallel from input_fasta except where noted):
 *   1. CheckV          — genome completeness + quality assessment
 *   2. geNomad         — virus classification + gene prediction
 *   3. RdRPCATCH       — RdRP detection (optional, default: true)
 *   4. BLASTn × N dbs  — nucleotide search (optional, default: true)
 *      └─ BLASTani       — pairwise ANI from BLASTn output
 *   5. DIAMOND         — protein search vs NCBI-nr (optional, default: false)
 *                        uses geNomad-predicted proteins
 *   6. RNAVirHost      — host prediction (optional, default: true)
 *                        requires CheckV + geNomad output
 *   7. MERGE_ANNOTATIONS — combine all outputs per vOTU
 *   8. GUESS_HOST      — ensemble LLM host prediction (optional, default: false)
 *   9. ANNOTATION_SUMMARY — per-step counts report
 *
 * BLASTn databases are specified via a CSV file (--blastn_dbs):
 *   name,path
 *   IMGVR,/path/to/IMGVR5_UViG.fna
 *   CHVD,/path/to/CHVD_virus_sequences.fasta
 *
 * Usage:
 *   nextflow run main.nf -profile cluster,conda \
 *       --input_fasta   /path/to/votus.fasta \
 *       --outdir        /path/to/results \
 *       --checkvdb      /path/to/checkv-db-v1.5 \
 *       --genomad_db    /path/to/genomad_db \
 *       --blastn_dbs    /path/to/blastn_dbs.csv
 */

nextflow.enable.dsl = 2

include { CHECKV             } from './modules/checkv'
include { GENOMAD            } from './modules/genomad'
include { RDRPCATCH          } from './modules/rdrpcatch'
include { BLASTN             } from './modules/blastn_tools'
include { BLASTANI           } from './modules/blastn_tools'
include { DIAMOND            } from './modules/diamond'
include { RNAVIRHOST         } from './modules/rnavirhost'
include { PREPARE_ICTV       } from './modules/summary'
include { MERGE_ANNOTATIONS  } from './modules/summary'
include { GUESS_HOST         } from './modules/summary'
include { ANNOTATION_SUMMARY } from './modules/summary'

workflow {

    // -----------------------------------------------------------------------
    // PREPARE_ICTV mode — run with --prepare_ictv true to generate the ICTV
    // reference file without running the full annotation workflow.
    //
    // Usage (auto-download):
    //   nextflow run main.nf --prepare_ictv true --outdir /path/to/ref_data
    //
    // Usage (manual download):
    //   nextflow run main.nf --prepare_ictv true \
    //       --outdir /path/to/ref_data \
    //       --ictv_raw /path/to/ictv_VirusPropertiesByFamily.tsv
    // -----------------------------------------------------------------------
    if ( params.prepare_ictv ) {
        raw_tsv = params.ictv_raw
            ? Channel.fromPath(params.ictv_raw, checkIfExists: !workflow.stubRun)
            : Channel.of(file('NO_FILE_ICTV_RAW'))
        PREPARE_ICTV(raw_tsv)
        return
    }

    // -----------------------------------------------------------------------
    // Parameter validation
    // -----------------------------------------------------------------------
    def errors = []
    if (!params.input_fasta)   errors << "  --input_fasta is required"
    if (!params.outdir)        errors << "  --outdir is required"
    if (!params.checkvdb)      errors << "  --checkvdb is required"
    if (!params.genomad_db)    errors << "  --genomad_db is required"
    // Nextflow may pass boolean CLI params as strings; normalise to boolean
    def run_rdrpcatch  = params.run_rdrpcatch  instanceof Boolean ? params.run_rdrpcatch  : params.run_rdrpcatch.toString()  != 'false'
    def run_blastn     = params.run_blastn     instanceof Boolean ? params.run_blastn     : params.run_blastn.toString()     != 'false'
    def run_diamond    = params.run_diamond    instanceof Boolean ? params.run_diamond    : params.run_diamond.toString()    != 'false'
    def run_rnavirhost = params.run_rnavirhost instanceof Boolean ? params.run_rnavirhost : params.run_rnavirhost.toString() != 'false'
    def run_guess_host = params.run_guess_host instanceof Boolean ? params.run_guess_host : params.run_guess_host.toString() != 'false'

    if (run_rdrpcatch  && !params.rdrpcatch_db)  errors << "  --rdrpcatch_db is required when --run_rdrpcatch=true"
    if (run_blastn     && !params.blastn_dbs)    errors << "  --blastn_dbs is required when --run_blastn=true"
    if (run_diamond    && !params.diamond_db)    errors << "  --diamond_db is required when --run_diamond=true"
    if (errors) {
        log.error "Parameter errors:\n" + errors.join("\n")
        System.exit(1)
    }

    // Validate paths exist (skip during stub runs)
    if (!workflow.stubRun) {
        [params.input_fasta, params.checkvdb, params.genomad_db].each { p ->
            if (p && !file(p).exists()) error "Path not found: ${p}"
        }
        if (run_rdrpcatch && !file(params.rdrpcatch_db).exists())
            error "rdrpcatch_db not found: ${params.rdrpcatch_db}"
        if (run_blastn && !file(params.blastn_dbs).exists())
            error "blastn_dbs CSV not found: ${params.blastn_dbs}"
        if (run_diamond && !file(params.diamond_db).exists())
            error "diamond_db not found: ${params.diamond_db}"
    }

    log.info """
    ============================================
     wvdb_annotate: viral genome annotation
    ============================================
     input_fasta        : ${params.input_fasta}
     outdir             : ${params.outdir}
     run_rdrpcatch      : ${params.run_rdrpcatch}
     run_blastn         : ${params.run_blastn}
     run_diamond        : ${params.run_diamond}
     run_rnavirhost     : ${params.run_rnavirhost}
     run_guess_host     : ${params.run_guess_host}
    ============================================
    """.stripIndent()

    def check      = !workflow.stubRun
    input_fasta    = file(params.input_fasta, checkIfExists: check)
    checkvdb       = file(params.checkvdb,    checkIfExists: check)
    genomad_db     = file(params.genomad_db,  checkIfExists: check)
    diamond_db     = params.diamond_db ?: ''  // val — not staged

    // -----------------------------------------------------------------------
    // Step 1 — CheckV
    // -----------------------------------------------------------------------
    CHECKV(
        input_fasta,
        checkvdb
    )

    // -----------------------------------------------------------------------
    // Step 2 — geNomad
    // -----------------------------------------------------------------------
    GENOMAD(
        input_fasta,
        genomad_db
    )

    // -----------------------------------------------------------------------
    // Step 3 — RdRPCATCH (optional)
    // -----------------------------------------------------------------------
    if ( run_rdrpcatch ) {
        rdrpcatch_db = file(params.rdrpcatch_db, checkIfExists: check)
        RDRPCATCH(
            input_fasta,
            rdrpcatch_db
        )
        rdrpcatch_tsv = RDRPCATCH.out.results_tsv.flatten().first()
    } else {
        rdrpcatch_tsv = Channel.of(file('NO_FILE_RDRPCATCH'))
    }

    // -----------------------------------------------------------------------
    // Step 4 — BLASTn × N databases (optional, parallel)
    // blastn_dbs CSV format: name,path  (one database per row, header required)
    // -----------------------------------------------------------------------
    if ( run_blastn ) {
        // Parse CSV → channel of [db_name, db_path] per row
        // BLAST db paths kept as strings (val) to avoid staging index file issues
        blastn_dbs_ch = Channel
            .fromPath(params.blastn_dbs)
            .splitCsv(header: true)
            .map { row -> tuple(row.name, row.path) }

        BLASTN(
            input_fasta,
            blastn_dbs_ch.map { name, path -> name },
            blastn_dbs_ch.map { name, path -> path }
        )

        BLASTANI(
            BLASTN.out.blastn_result
        )

        // Collect all ANI TSVs for MERGE_ANNOTATIONS
        // merge_annotations.py receives them as a flat list staged in the work dir
        blastn_ani_tsvs = BLASTANI.out.ani_result
            .map { db_name, tsv -> tsv }
            .collect()
    } else {
        blastn_ani_tsvs = Channel.of(file('NO_FILE_BLASTN'))
    }

    // -----------------------------------------------------------------------
    // Step 5 — DIAMOND protein search (optional)
    // Uses predicted proteins from geNomad
    // -----------------------------------------------------------------------
    if ( run_diamond ) {
        DIAMOND(
            GENOMAD.out.proteins_faa,
            diamond_db
        )
        diamond_tsv = DIAMOND.out.diamond_tsv
    } else {
        diamond_tsv = Channel.of(file('NO_FILE_DIAMOND'))
    }

    // -----------------------------------------------------------------------
    // Step 6 — RNAVirHost host prediction (optional)
    // Sequential: requires CheckV + geNomad output
    // -----------------------------------------------------------------------
    if ( run_rnavirhost ) {
        RNAVIRHOST(
            input_fasta,
            CHECKV.out.quality_summary,
            GENOMAD.out.virus_summary
        )
        rnavirhost_tsv = RNAVIRHOST.out.results_tsv.flatten().first()
    } else {
        rnavirhost_tsv = Channel.of(file('NO_FILE_RNAVIRHOST'))
    }

    // -----------------------------------------------------------------------
    // Step 7 — Merge all annotations
    // -----------------------------------------------------------------------
    MERGE_ANNOTATIONS(
        CHECKV.out.quality_summary,
        GENOMAD.out.virus_summary,
        file(params.ictv_fam),
        blastn_ani_tsvs,
        rdrpcatch_tsv,
        rnavirhost_tsv
    )

    // -----------------------------------------------------------------------
    // Step 8 — Guess host (optional, sequential after merge)
    // -----------------------------------------------------------------------
    if ( run_guess_host ) {
        GUESS_HOST(
            MERGE_ANNOTATIONS.out.merged_tsv,
            input_fasta
        )
        host_tsv = GUESS_HOST.out.host_tsv
    } else {
        host_tsv = Channel.of(file('NO_FILE_GUESS_HOST'))
    }

    // -----------------------------------------------------------------------
    // Step 9 — Annotation summary report
    // -----------------------------------------------------------------------
    ANNOTATION_SUMMARY(
        input_fasta,
        CHECKV.out.quality_summary,
        GENOMAD.out.virus_summary,
        MERGE_ANNOTATIONS.out.merged_tsv,
        host_tsv
    )
}

