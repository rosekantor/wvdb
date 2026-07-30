/*
 * modules/collect.nf
 * Step 1: Collect filtered virus and provirus FASTAs from a list of paths,
 * merge into a single filtered_all.fasta, and compute per-sequence lengths.
 *
 * Input: a plain text file with one absolute FASTA path per line.
 * Lines containing "provirus" are treated as provirus FASTAs.
 * Lines containing "virus" but not "provirus" are treated as virus FASTAs.
 * Empty lines and lines starting with '#' are ignored.
 *
 * Example fasta_list.txt:
 *   /p/vast1/.../batch1/assembly/sample_all/filtered_virus_all.fasta
 *   /p/vast1/.../batch1/assembly/sample_all/filtered_provirus_all.fasta
 *   /p/vast1/.../batch2/assembly/sample_all/filtered_virus_all.fasta
 *   /p/vast1/.../batch2/assembly/sample_all/filtered_provirus_all.fasta
 */

process COLLECT_GENOMES {
    label 'cpu_low'

    publishDir "${params.outdir}/1_collect", mode: 'copy', enabled: !workflow.stubRun

    input:
    path fasta_list   // text file of absolute FASTA paths, one per line

    output:
    path "filtered_all.fasta",          emit: merged_fasta
    path "filtered_all.length.txt",     emit: lengths
    path "filtered_virus_all.fasta",    emit: virus_fasta
    path "filtered_provirus_all.fasta", emit: provirus_fasta

    script:
    """
    virus_files=()
    provirus_files=()

    while IFS= read -r line || [[ -n "\$line" ]]; do
        [[ -z "\$line" || "\$line" == \\#* ]] && continue
        if [[ "\$line" == *provirus* ]]; then
            provirus_files+=( "\$line" )
        elif [[ "\$line" == *virus* ]]; then
            virus_files+=( "\$line" )
        else
            echo "Warning: could not classify as virus or provirus: \$line" >&2
        fi
    done < "${fasta_list}"

    echo "Found \${#virus_files[@]} virus FASTA(s) and \${#provirus_files[@]} provirus FASTA(s)" >&2

    [[ \${#virus_files[@]} -gt 0 || \${#provirus_files[@]} -gt 0 ]] || {
        echo "ERROR: no valid FASTA paths found in ${fasta_list}" >&2
        exit 1
    }

    if [[ \${#virus_files[@]} -gt 0 ]]; then
        cat "\${virus_files[@]}" > filtered_virus_all.fasta
    else
        echo "Warning: no virus FASTAs — filtered_virus_all.fasta will be empty" >&2
        touch filtered_virus_all.fasta
    fi

    if [[ \${#provirus_files[@]} -gt 0 ]]; then
        cat "\${provirus_files[@]}" > filtered_provirus_all.fasta
    else
        echo "Warning: no provirus FASTAs — filtered_provirus_all.fasta will be empty" >&2
        touch filtered_provirus_all.fasta
    fi

    cat filtered_virus_all.fasta filtered_provirus_all.fasta > filtered_all.fasta

    [[ -s filtered_all.fasta ]] || {
        echo "ERROR: merged filtered_all.fasta is empty" >&2
        exit 1
    }

    seqkit fx2tab -nl filtered_all.fasta > filtered_all.length.txt

    echo "Total sequences: \$(grep -c '^>' filtered_all.fasta)" >&2
    """

    stub:
    """
    touch filtered_all.fasta filtered_all.length.txt filtered_virus_all.fasta filtered_provirus_all.fasta
    """
}
