/*
 * modules/protein_tools.nf
 * Processes for protein characterization.
 *
 * Process list:
 *   DIAMOND      — protein similarity search against a reference database
 *   HMMSEARCH    — HMM profile search against predicted proteins
 *   PROTEIN_SUMMARY — per-contig protein hit counts for MERGE_ANNOTATIONS
 *
 * Both DIAMOND and HMMSEARCH are called once per database/profile via
 * channels derived from diamond_dbs.csv and hmm_profiles.csv respectively.
 *
 * DIAMOND databases are passed as val (not path) to prevent staging index
 * file issues (same pattern as BLASTn databases).
 *
 * HMM profile files are passed as path since they are single files without
 * associated index files.
 */

// ---------------------------------------------------------------------------
// DIAMOND
// Protein similarity search against a prebuilt DIAMOND database.
// db_name scopes output filenames; db_path passed as val to avoid staging.
//
// Output format: outfmt 6 + stitle (+ staxids if params.diamond_taxonmap=true)
// Default filter: bitscore >= params.diamond_min_bitscore (default 50)
// ---------------------------------------------------------------------------
process DIAMOND {
    label 'cpu_high'

    tag "${db_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${db_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path  proteins_faa            // geNomad predicted proteins .faa
    val   db_name                 // short label from diamond_dbs.csv
    val   db_path                 // path to .dmnd database (val not path)

    output:
    tuple val(db_name), path("${db_name}.diamond.tsv"), emit: diamond_result

    script:
    def taxids_field = params.diamond_taxonmap ? " staxids" : ""
    def outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle${taxids_field}"
    """
    diamond blastp \\
        --query        "${proteins_faa}" \\
        --db           "${db_path}" \\
        --out          "${db_name}.diamond.tsv" \\
        --outfmt       ${outfmt} \\
        --evalue       ${params.diamond_evalue} \\
        --min-score    ${params.diamond_min_bitscore} \\
        --max-target-seqs ${params.diamond_max_targets} \\
        --threads      ${task.cpus} \\
        --more-sensitive
    """

    stub:
    """
    touch ${db_name}.diamond.tsv
    """
}

// ---------------------------------------------------------------------------
// HMMSEARCH
// HMM profile search against predicted proteins.
// Uses --cut_tc (trusted cutoff) if available in the profile; falls back to
// --domE params.hmm_evalue if --cut_tc causes an error (profiles without TC).
//
// Output: domain table (--domtblout) parsed by parse_hmmsearch.py,
// plus the raw --domtblout for archival.
// ---------------------------------------------------------------------------
process HMMSEARCH {
    label 'cpu_high'

    tag "${profile_name}"

    publishDir { "${params.outdir}/${task.ext.publish_dir}/${profile_name}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path  proteins_faa            // geNomad predicted proteins .faa
    val   profile_name            // short label from hmm_profiles.csv
    path  hmm_profile             // HMM profile file (.hmm)

    output:
    tuple val(profile_name), path("${profile_name}.hmmsearch.tsv"),     emit: hmmsearch_result
    tuple val(profile_name), path("${profile_name}.hmmsearch.domtbl"),  emit: domtbl_raw

    script:
    """
    # Try --cut_tc first (trusted cutoffs embedded in profile)
    # Fall back to E-value filter if profile has no TC thresholds
    if hmmsearch \\
            --domtblout "${profile_name}.hmmsearch.domtbl" \\
            --cut_tc \\
            --cpu ${task.cpus} \\
            --noali \\
            "${hmm_profile}" \\
            "${proteins_faa}" > /dev/null 2> hmm_err.txt; then
        echo "hmmsearch completed with --cut_tc" >&2
    elif grep -qi "no TC cutoff" hmm_err.txt || grep -qi "thresholds unavailable" hmm_err.txt; then
        echo "No TC cutoffs in profile — falling back to E-value filter (${params.hmm_evalue})" >&2
        hmmsearch \\
            --domtblout "${profile_name}.hmmsearch.domtbl" \\
            --domE ${params.hmm_evalue} \\
            --cpu ${task.cpus} \\
            --noali \\
            "${hmm_profile}" \\
            "${proteins_faa}" > /dev/null
    else
        echo "hmmsearch failed:" >&2
        cat hmm_err.txt >&2
        exit 1
    fi

    parse_hmmsearch.py \\
        --domtbl "${profile_name}.hmmsearch.domtbl" \\
        --out    "${profile_name}.hmmsearch.tsv"
    """

    stub:
    """
    touch ${profile_name}.hmmsearch.tsv ${profile_name}.hmmsearch.domtbl
    """
}

// ---------------------------------------------------------------------------
// PROTEIN_SUMMARY
// Per-contig protein hit counts across all diamond and hmmsearch results.
// Adds n_proteins_total and n_proteins_with_hit_<db> columns for each db.
// ---------------------------------------------------------------------------
process PROTEIN_SUMMARY {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path  proteins_faa            // geNomad predicted proteins .faa
    path  diamond_tsvs            // collected diamond result TSVs (or NO_FILE_DIAMOND)
    path  hmmsearch_tsvs          // collected hmmsearch result TSVs (or NO_FILE_HMMSEARCH)

    output:
    path "protein_summary.tsv", emit: protein_summary_tsv

    script:
    def diamond_arg   = !diamond_tsvs.every  { it.name.startsWith('NO_FILE') } ? "--diamond-dir   ." : ""
    def hmmsearch_arg = !hmmsearch_tsvs.every { it.name.startsWith('NO_FILE') } ? "--hmmsearch-dir ." : ""
    """
    summarize_protein_hits.py \\
        --faa     "${proteins_faa}" \\
        ${diamond_arg} \\
        ${hmmsearch_arg} \\
        --out     protein_summary.tsv
    """

    stub:
    """
    touch protein_summary.tsv
    """
}
