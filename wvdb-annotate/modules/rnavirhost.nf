/*
 * modules/rnavirhost.nf
 * RNAVIRHOST — host prediction for RNA viruses.
 *
 * Builds a consensus viral order classification from RdRPCATCH + geNomad
 * outputs (RdRPCATCH preferred; geNomad fallback; 'Unclassified' if neither),
 * then runs `rnavirhost predict`.
 *
 * RNAVirHost is available on bioconda:
 *   conda install -c bioconda rnavirhost
 * It is included in the main wvdb-annotate conda environment.
 *
 * Note: RNAVirHost is designed for RNA viruses. Sequences without protein-
 * coding genes (detected by Prodigal) fall back to BLASTn-based host
 * assignment. DNA viruses will mostly receive 'unclassified' evidence.
 *
 * The --force flag is passed to run_rnavirhost.py so Nextflow -resume
 * does not fail if the output directory already exists in the work dir.
 */

process RNAVIRHOST {
    label 'cpu_medium'

    publishDir { "${params.outdir}/${task.ext.publish_dir}" }, mode: 'copy',
        enabled: !workflow.stubRun

    input:
    path input_fasta
    path checkv_quality           // checkv_out/quality_summary.tsv
    path genomad_virus_summary    // genomad <prefix>_virus_summary.tsv
    path rdrpcatch_tsv            // rdrpcatch annotated output TSV

    output:
    path "rnavirhost_out/",                             emit: rnavirhost_dir
    path "rnavirhost_out/predict/result.csv",           emit: result_csv,  optional: true
    path "rnavirhost_out/result.csv",                   emit: result_csv2, optional: true
    path "rnavirhost_consensus_orders.csv",             emit: orders_csv

    script:
    """
    run_rnavirhost.py \\
        --genomad  "${genomad_virus_summary}" \\
        --checkv   "${checkv_quality}" \\
        --rdrpcatch "${rdrpcatch_tsv}" \\
        -f "${input_fasta}" \\
        -O rnavirhost_consensus_orders.csv \\
        -o rnavirhost_out \\
        --force
    """

    stub:
    """
    mkdir -p rnavirhost_out/predict
    touch rnavirhost_out/predict/result.csv
    touch rnavirhost_out/result.csv
    touch rnavirhost_consensus_orders.csv
    """
}
