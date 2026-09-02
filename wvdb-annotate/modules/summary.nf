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
// ---------------------------------------------------------------------------
// CATEGORIZE_HOST_TERMS
// LLM-based categorization of free-text host/isolation-source strings from
// core_nt BLAST hits. Separated from GUESS_HOST so the slow, API-cost-
// incurring categorization step is cached and only reruns when the merged
// annotation table's unique terms actually change — the decision-tree logic
// downstream can be freely rerun/tuned without re-querying an LLM.
//
// .env file is read directly from the filesystem by categorize_host_terms.py
// (params.llm_env_file, passed as a plain string, not staged as a Nextflow
// path input) — this keeps the API key out of Nextflow's work/trace dirs.
// ---------------------------------------------------------------------------
process CATEGORIZE_HOST_TERMS {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path merged_annotations       // output of MERGE_ANNOTATIONS
    path host_dict_cache          // previous host_dict.tsv (or NO_FILE_HOST_DICT)
    path isolation_dict_cache     // previous isolation_dict.tsv (or NO_FILE_ISO_DICT)

    output:
    path "host_dict.tsv",      emit: host_dict
    path "isolation_dict.tsv", emit: isolation_dict
    path "rejected_terms.tsv", emit: rejected_terms

    script:
    def host_cache_arg = !host_dict_cache.name.startsWith('NO_FILE') \
        ? "--host-dict ${host_dict_cache}" : ""
    def iso_cache_arg  = !isolation_dict_cache.name.startsWith('NO_FILE') \
        ? "--isolation-dict ${isolation_dict_cache}" : ""
    def base_url_arg   = params.llm_base_url \
        ? "--llm-base-url ${params.llm_base_url}" : ""
    def api_key_env_arg = params.llm_api_key_env_var \
        ? "--llm-api-key-env-var ${params.llm_api_key_env_var}" : ""
    """
    categorize_host_terms.py \\
        --merged              "${merged_annotations}" \\
        --env-file             "${params.llm_env_file}" \\
        --provider              ${params.llm_provider} \\
        --model                 ${params.llm_model} \\
        ${base_url_arg} \\
        ${api_key_env_arg} \\
        ${host_cache_arg} \\
        ${iso_cache_arg} \\
        --out-host-dict         host_dict.tsv \\
        --out-isolation-dict    isolation_dict.tsv \\
        --out-rejected-terms    rejected_terms.tsv
    """

    stub:
    """
    touch host_dict.tsv isolation_dict.tsv rejected_terms.tsv
    """
}

// ---------------------------------------------------------------------------
// GUESS_HOST
// Ensemble host determination from ICTV, RNAVirHost, and (if categorized)
// core_nt BLAST hit metadata. Pure pandas decision-tree logic, no LLM
// dependency — safe to rerun freely.
// Optional: controlled by params.run_guess_host (default: false).
// ---------------------------------------------------------------------------
process GUESS_HOST {
    label 'cpu_low'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path merged_annotations    // output of MERGE_ANNOTATIONS
    path host_dict             // from CATEGORIZE_HOST_TERMS (or NO_FILE)
    path isolation_dict        // from CATEGORIZE_HOST_TERMS (or NO_FILE)

    output:
    path "final_annotations.tsv", emit: host_tsv
    path "manual_review.tsv",    emit: manual_review_tsv

    script:
    def host_dict_arg = !host_dict.name.startsWith('NO_FILE') \
        ? "--host-dict ${host_dict}" : ""
    def iso_dict_arg   = !isolation_dict.name.startsWith('NO_FILE') \
        ? "--isolation-dict ${isolation_dict}" : ""
    """
    guess_host.py \\
        --merged              "${merged_annotations}" \\
        ${host_dict_arg} \\
        ${iso_dict_arg} \\
        --out                  final_annotations.tsv \\
        --out-manual-review    manual_review.tsv \\
        --other-threshold      ${params.guess_host_other_threshold}
    """

    stub:
    """
    touch final_annotations.tsv manual_review.tsv
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
    path host_predictions        // final_annotations.tsv (or NO_FILE)

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
