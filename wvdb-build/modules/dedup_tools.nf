/*
 * modules/dedup_tools.nf
 * DEDUPLICATE_GENOMES — detect and correct whole-genome duplications
 * (tandem or inverted) flagged by CheckV's kmer_freq metric.
 *
 * Runs right before the final reclustering step (step 5), on the merged
 * complete representatives from cluster_trim + blast_trim. This is a
 * common assembly artifact — an assembler occasionally outputs a genome
 * as two tandem or inverted-repeat copies concatenated together — and
 * CheckV's kmer_freq (~2.0) flags it directly.
 *
 * Detection is alignment-based (nucmer): a short probe from the start of
 * each flagged sequence is aligned against the full sequence. Both
 * orientations are considered since nucmer searches both strands by
 * default. Only sequences with clean, well-separated, high-identity,
 * evenly-spaced repeat copies are auto-corrected; anything ambiguous
 * (e.g. a rotated/irregular repeat) is left untouched and flagged for
 * manual review rather than risking an incorrect edit.
 *
 * See bin/detect_genome_duplication.py for the full algorithm and
 * bin/detect_genome_duplication.py's module docstring for validation
 * details (tested against synthetic tandem, inverted, rotated, and
 * unrelated-sequence cases before being used on real data).
 */

process DEDUPLICATE_GENOMES {
    label 'cpu_medium'

    publishDir "${params.outdir}/dedup", mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path input_fasta
    path checkv_quality

    output:
    path "deduplicated_all.fasta",   emit: deduplicated_fasta
    path "corrected_only.fasta",     emit: corrected_fasta
    path "flagged_for_review.fasta", emit: flagged_fasta
    path "dedup_report.tsv",         emit: report_tsv
    path "plots/",                   emit: plots_dir

    script:
    """
    detect_genome_duplication.py \\
        --fasta         "${input_fasta}" \\
        --checkv        "${checkv_quality}" \\
        --outdir        . \\
        --threshold     ${params.dedup_kmer_threshold} \\
        --min-identity  ${params.dedup_min_identity} \\
        --min-coverage  ${params.dedup_min_coverage}
    """

    stub:
    """
    mkdir -p plots
    touch deduplicated_all.fasta corrected_only.fasta flagged_for_review.fasta dedup_report.tsv
    """
}
