/*
 * modules/summary.nf
 * Annotation merging, host prediction, and summary report processes.
 *
 * Process list:
 *   MERGE_ANNOTATIONS  — combine all annotation outputs into one TSV per vOTU
 *   GUESS_HOST         — ensemble host prediction (optional LLM step)
 *   ANNOTATION_SUMMARY — per-step counts report (markdown + TSV)
 *
 * MERGE_ANNOTATIONS and GUESS_HOST are stubbed pending upload of
 * merge_annotations.py and guess_host.py. Interfaces designed to match
 * expected script signatures; update when scripts are uploaded.
 *
 * All optional inputs use NO_FILE sentinel so the process runs correctly
 * whether or not upstream optional steps were enabled.
 */

// ---------------------------------------------------------------------------
// MERGE_ANNOTATIONS
// Combines outputs from CheckV, geNomad, RdRPCATCH, BLASTn, DIAMOND,
// and RNAVirHost into a single per-vOTU annotation TSV.
// ---------------------------------------------------------------------------
process MERGE_ANNOTATIONS {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path checkv_quality           // checkv_out/quality_summary.tsv
    path genomad_virus_summary    // genomad virus_summary.tsv
    path ictv_fam                 // ref_data/ictv_families.tsv
    path blastn_ani_tsvs          // <db_name>.blastn.ani.tsv files (collected, or NO_FILE_BLASTN)
    path rdrpcatch_tsv            // rdrpcatch annotated TSV (or NO_FILE_RDRPCATCH)
    path rnavirhost_tsv           // rnavirhost result.csv (or NO_FILE_RNAVIRHOST)

    output:
    path "merged_annotations.tsv",  emit: merged_tsv
    path "host_lineage_cache.tsv",  emit: host_lineage_cache, optional: true

    script:
    def blastn_arg  = !(blastn_ani_tsvs instanceof List ? blastn_ani_tsvs[0] : blastn_ani_tsvs).name.startsWith('NO_FILE') ? "--blastn-dir ." : ""
    def rdrp_arg    = !rdrpcatch_tsv.name.startsWith('NO_FILE')  ? "--rdrpcatch  ${rdrpcatch_tsv}"  : ""
    def rnavirh_arg = !rnavirhost_tsv.name.startsWith('NO_FILE') ? "--rnavirhost ${rnavirhost_tsv}" : ""
    """
    merge_annotations.py \\
        --checkv       "${checkv_quality}" \\
        --genomad      "${genomad_virus_summary}" \\
        --ictv-fam     "${ictv_fam}" \\
        --entrez-email "${params.entrez_email}" \\
        ${blastn_arg} \\
        ${rdrp_arg} \\
        ${rnavirh_arg} \\
        --host-lineage-cache host_lineage_cache.tsv \\
        --out          merged_annotations.tsv
    """

    stub:
    """
    touch merged_annotations.tsv host_lineage_cache.tsv
    """
}

// ---------------------------------------------------------------------------
// GUESS_HOST
// Ensemble host prediction combining BLASTn taxonomy, geNomad predictions,
// RNAVirHost output, and optionally an LLM API call.
// Optional: controlled by params.run_guess_host (default: false).
// ---------------------------------------------------------------------------
process GUESS_HOST {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path merged_annotations    // output of MERGE_ANNOTATIONS
    path input_fasta

    output:
    path "host_predictions.tsv", emit: host_tsv

    script:
    """
    guess_host.py \\
        --annotations "${merged_annotations}" \\
        --fasta       "${input_fasta}" \\
        --out         host_predictions.tsv
    """

    stub:
    """
    touch host_predictions.tsv
    """
}

// ---------------------------------------------------------------------------
// ANNOTATION_SUMMARY
// Counts sequences annotated at each step and writes a report.
// ---------------------------------------------------------------------------
process ANNOTATION_SUMMARY {
    label 'cpu_low'

    publishDir "${params.outdir}", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path input_fasta
    path checkv_quality
    path genomad_virus_summary
    path merged_annotations
    path host_predictions        // host_predictions.tsv (or NO_FILE)

    output:
    path "annotation_summary.tsv", emit: summary_tsv
    path "annotation_summary.md",  emit: summary_md

    script:
    def host_arg = !host_predictions.name.startsWith('NO_FILE') \
        ? "--host-predictions ${host_predictions}" : ""
    """
    annotation_summary.py \\
        --fasta      "${input_fasta}" \\
        --checkv     "${checkv_quality}" \\
        --genomad    "${genomad_virus_summary}" \\
        --merged     "${merged_annotations}" \\
        ${host_arg} \\
        --out-tsv    annotation_summary.tsv \\
        --out-md     annotation_summary.md \\
        --run-date   "\$(date +%Y-%m-%d)"
    """

    stub:
    """
    touch annotation_summary.tsv annotation_summary.md
    """
}

// ---------------------------------------------------------------------------
// PREPARE_ICTV
// One-time utility process to process a manually saved ICTV family HTML page.
// To obtain the input:
//   1. Open https://ictv.global/virus-properties in your browser
//   2. Set 'Items per page' to 'All'
//   3. File → Save Page As → save as HTML
//   4. Pass the saved file via --ictv_raw
// Run: nextflow run main.nf --prepare_ictv true \
//          --ictv_raw /path/to/saved.html --outdir /path/to/ref_data
// Re-run when a new ICTV release is available.
// ---------------------------------------------------------------------------
process PREPARE_ICTV {
    label 'cpu_low'

    publishDir "${params.outdir}/ref_data", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path ictv_raw_tsv   // manually downloaded TSV, or NO_FILE to auto-download

    output:
    path "ictv_families.tsv", emit: ictv_tsv

    script:
    def in_arg = ictv_raw_tsv.name.startsWith('NO_FILE') ? '' : "--in ${ictv_raw_tsv}"
    """
    prepare_ictv.py ${in_arg} --out ictv_families.tsv
    """

    stub:
    """
    touch ictv_families.tsv
    """
}
