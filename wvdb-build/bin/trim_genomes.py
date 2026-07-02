#!/usr/bin/env python3
"""
trim_genomes.py — trim query sequences using nucmer pairwise alignments.

Cluster mode  (-c): aligns rank1 vs rank2/3 pairs from parse_clusters.py.
BLAST mode    (-a): aligns query vs best BLAST hit from blastani_nayfach.py.

Key improvement over the original version:
  In cluster mode, the input FASTA (-f) is already the small pre-extracted
  trimming_candidates.fasta (rank1/2/3 sequences only, produced by
  EXTRACT_TRIMMING_SEQS in the Nextflow pipeline). seqkit builds its faidx
  index on this small file once, making all per-pair extractions fast.
  In blast mode, the query FASTA is trimming_candidates.fasta and the
  target FASTA is the BLAST database FASTA, unchanged from before.
"""

import argparse
import concurrent.futures
from glob import glob
import os
import pandas as pd
import re
import shutil
import subprocess
import sys


# ---------------------------------------------------------------------------
# Initialisation
# ---------------------------------------------------------------------------

def initialize(query_fasta, outdir):
    """Check for dependencies and required files. Make output subdirectories."""
    REQUIRED_PROGS = ["seqkit", "nucmer", "show-coords"]
    missing = [p for p in REQUIRED_PROGS if shutil.which(p) is None]
    if missing:
        print(
            "Error: required program(s) not found in PATH: {}\n"
            "Please install them or add them to your PATH.\n".format(", ".join(missing)),
            file=sys.stderr
        )
        sys.exit(1)

    if not os.path.isfile(query_fasta):
        print(f"Error: fasta file not found: {query_fasta}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, 'pairs'), exist_ok=True)
    os.makedirs(os.path.join(outdir, 'aln'),   exist_ok=True)


def build_faidx(fasta_path, threads):
    """
    Build a seqkit .fai index for fasta_path if one does not already exist.
    Skipping the rebuild saves significant time for large databases (100s of GB).
    Returns True on success, False on failure.
    """
    fai_path = fasta_path + ".fai"
    if os.path.isfile(fai_path):
        print(f"  faidx already exists, skipping: {fai_path}", file=sys.stderr)
        return True
    cmd = ["seqkit", "faidx", fasta_path, "--threads", str(threads)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error building faidx for {fasta_path}: {result.stderr}", file=sys.stderr)
        return False
    return True


# ---------------------------------------------------------------------------
# Pre-extract all target sequences in one seqkit call
# ---------------------------------------------------------------------------

def extract_targets(pairs_df, target_fasta, fastadir, threads):
    """
    Extract all unique target sequences from target_fasta in a single
    seqkit faidx call, writing them to a small per-targets FASTA.
    This avoids N individual subprocess calls against the large database
    during parallel nucmer alignment.
    Returns path to the extracted targets FASTA, or None on failure.
    """
    targets = pairs_df['tname'].unique().tolist()
    targets_fasta = os.path.join(fastadir, 'all_targets.fasta')

    # Write target IDs to a temp file for seqkit grep
    ids_file = os.path.join(fastadir, 'target_ids.txt')
    with open(ids_file, 'w') as f:
        for t in targets:
            f.write(t + '\n')

    cmd = ['seqkit', 'grep', '-f', ids_file, target_fasta,
           '-o', targets_fasta, '--threads', str(threads)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error extracting target sequences: {result.stderr}", file=sys.stderr)
        return None

    print(f"  Pre-extracted {len(targets)} target sequences → {targets_fasta}",
          file=sys.stderr)

    # Build faidx on the small targets FASTA for fast per-pair access
    build_faidx(targets_fasta, threads)
    return targets_fasta


# ---------------------------------------------------------------------------
# BLAST mode pair extraction (unchanged from original)
# ---------------------------------------------------------------------------

def get_alns(blastani_file, qcov=85, tcov=85):
    """
    Load BLAST ANI results and derive query/target pairs DataFrame.
    Filters on qcov OR tcov >= threshold, keeps best hit per query.
    """
    ani_df = pd.read_csv(blastani_file, sep='\t')
    ani_df = ani_df[(ani_df.qcov >= qcov) | (ani_df.tcov >= tcov)]
    if ani_df.empty:
        print("Warning: no BLAST ANI hits passed coverage filters.", file=sys.stderr)
        return pd.DataFrame(columns=["qname", "tname"])
    idx = ani_df.groupby('qname')['qcov'].idxmax()
    return ani_df.loc[idx][['qname', 'tname']]


# ---------------------------------------------------------------------------
# Per-pair alignment (parallel worker)
# ---------------------------------------------------------------------------

def pull_aln(row_id, query, target, query_fasta, target_fasta, fastadir, alndir):
    """
    Extract one query and one target sequence, align with nucmer, get coords.
    The faidx index on query_fasta and target_fasta must already exist
    (built by build_faidx before the parallel executor starts).

    show-coords flags (unchanged from original):
      -r  sort by reference
      -c  include coverage
      -l  include sequence lengths
      -L 1000  minimum alignment length 1000 nt
      -I 90    minimum identity 90%
      -T  tab-delimited output
    """
    target_newname = re.split(r'\||\ |,', target)[0]
    qsuffix = f"{query}_{row_id}"
    tsuffix = f"{target_newname}_{row_id}"

    cmd = (
        f'seqkit faidx "{query_fasta}"  "{query}"  > "{fastadir}/{qsuffix}.query.fasta";  '
        f'seqkit faidx "{target_fasta}" "{target}" > "{fastadir}/{tsuffix}.target.fasta"; '
        f'nucmer -p "{alndir}/query_{qsuffix}" '
        f'    "{fastadir}/{qsuffix}.query.fasta" '
        f'    "{fastadir}/{tsuffix}.target.fasta"; '
        f'show-coords -r -c -l -L 1000 -I 90 -T '
        f'    "{alndir}/query_{qsuffix}.delta" '
        f'    > "{alndir}/query_{qsuffix}.coords"'
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error for {query}: {result.stderr}", file=sys.stderr)
    return result.returncode


# ---------------------------------------------------------------------------
# Parse nucmer output → BED (unchanged from original)
# ---------------------------------------------------------------------------

def parse_nucmer(alndir):
    """
    Read all nucmer .coords files in alndir and produce a BED DataFrame.
    Takes the alignment with highest reference coverage per query sequence.
    """
    files = glob(f'{alndir}/*.coords')

    if not files:
        print(f"No .coords files found in {alndir}, cannot create BED.", file=sys.stderr)
        return pd.DataFrame(columns=['ref_name', 'start1', 'end1'])

    headers = [
        'start1', 'end1', 'start2', 'end2', 'len1', 'len2', 'pid',
        'len_ref', 'len_query', 'cov_ref', 'cov_query', 'ref_name', 'query_name'
    ]
    dfs = []
    for file in files:
        df = pd.read_csv(file, sep='\t', names=headers, skiprows=4)
        if not df.empty:
            dfs.append(df)

    if not dfs:
        print(f"Warning: all coords files in {alndir} are empty.", file=sys.stderr)
        return pd.DataFrame(columns=['ref_name', 'start1', 'end1'])

    aln_df = pd.concat(dfs, ignore_index=True)

    # Ensure start < end
    aln_df['start1'], aln_df['end1'] = (
        aln_df[['start1', 'end1']].min(axis=1),
        aln_df[['start1', 'end1']].max(axis=1)
    )

    # Keep alignment with highest reference coverage per query
    idx = aln_df.groupby('ref_name')['cov_ref'].idxmax()
    bed_df = aln_df.loc[idx][['ref_name', 'start1', 'end1']].copy()
    return bed_df


# ---------------------------------------------------------------------------
# Run seqkit subseq + clean headers (unchanged from original)
# ---------------------------------------------------------------------------

def run_trimming(query_fasta, outdir, threads):
    """
    Run seqkit subseq with the BED file to trim sequences, then
    clean up coordinate suffixes added by seqkit to FASTA headers.
    """
    bed_file    = os.path.join(outdir, "trimming.bed")
    temp_fasta  = os.path.join(outdir, "trimmed_temp.fasta")
    final_fasta = os.path.join(outdir, "trimmed.fasta")

    if not os.path.isfile(bed_file):
        print(f"Error: BED file '{bed_file}' not found.", file=sys.stderr)
        sys.exit(1)

    cmd_subseq = ["seqkit", "subseq", "--bed", bed_file, query_fasta,
                  "-j", str(threads)]
    try:
        with open(temp_fasta, "w") as fout:
            result = subprocess.run(cmd_subseq, stdout=fout,
                                    stderr=subprocess.PIPE, text=True)
    except OSError as e:
        print(f"Error writing temporary trimmed FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(f"Error running seqkit subseq: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # Remove coordinate suffixes added by seqkit (_start-end:.)
    pattern = re.compile(r'_\d+-\d+:\.')
    try:
        with open(temp_fasta) as fin, open(final_fasta, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    line = pattern.sub("", line)
                fout.write(line)
    except OSError as e:
        print(f"Error rewriting trimmed FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        os.remove(temp_fasta)
    except OSError:
        print(f"Warning: could not remove temp file '{temp_fasta}'", file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Trim query sequences using nucmer pairwise alignments."
    )
    parser.add_argument("-f", "--fasta", required=True,
                        help="Input FASTA (trimming_candidates.fasta in cluster mode, "
                             "or query FASTA in BLAST mode)")
    parser.add_argument("-o", "--outdir", required=True,
                        help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of parallel nucmer workers (default: 1)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-c", "--cluster_pairs",
                       help="Cluster pairs TSV from parse_clusters.py (cluster mode)")
    group.add_argument("-a", "--blast_ani",
                       help="BLAST ANI TSV from blastani_nayfach.py (BLAST mode)")

    parser.add_argument("-d", "--blast_db_fasta",
                        help="Reference FASTA for BLAST database (required in BLAST mode)")

    args = parser.parse_args()

    if args.blast_ani is not None and args.blast_db_fasta is None:
        parser.error("--blast_db_fasta is required when --blast_ani is used")
    if args.cluster_pairs is not None and args.blast_db_fasta is not None:
        parser.error("--blast_db_fasta is only used with --blast_ani")

    return args


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    threads     = args.threads
    query_fasta = args.fasta
    outdir      = args.outdir

    initialize(query_fasta, outdir)

    # --- Step 1: load pairs ---
    if args.cluster_pairs is not None:
        # Cluster mode: query_fasta is trimming_candidates.fasta
        # target_fasta is the same file (rank2/3 are in the same candidates FASTA)
        target_fasta = query_fasta
        if not os.path.isfile(args.cluster_pairs):
            print(f"Error: pairs file not found: {args.cluster_pairs}", file=sys.stderr)
            sys.exit(1)
        pairs_df = pd.read_csv(args.cluster_pairs, sep='\t', names=['qname', 'tname'])
    else:
        # BLAST mode: target_fasta is the external BLAST db FASTA
        target_fasta = args.blast_db_fasta
        for f in [args.blast_ani, target_fasta]:
            if not os.path.isfile(f):
                print(f"Error: file not found: {f}", file=sys.stderr)
                sys.exit(1)
        pairs_df = get_alns(args.blast_ani)

    if pairs_df.empty:
        print("No pairs to process — writing empty outputs.", file=sys.stderr)
        open(os.path.join(outdir, 'trimming.bed'),  'w').close()
        open(os.path.join(outdir, 'trimmed.fasta'), 'w').close()
        sys.exit(0)

    fastadir = os.path.join(outdir, 'pairs')
    alndir   = os.path.join(outdir, 'aln')

    # --- Step 2: prepare FASTAs for parallel per-pair extraction ---
    print(f"Building faidx index: {query_fasta}", file=sys.stderr)
    if not build_faidx(query_fasta, threads):
        sys.exit(1)

    if target_fasta != query_fasta:
        # BLAST mode: pre-extract all needed target sequences into a small
        # FASTA in one seqkit call. Workers then extract from this small file
        # rather than making individual calls to the large database.
        print(f"Pre-extracting target sequences from: {target_fasta}", file=sys.stderr)
        extracted = extract_targets(pairs_df, target_fasta, fastadir, threads)
        if extracted is None:
            sys.exit(1)
        target_fasta = extracted
    else:
        # Cluster mode: query and target are the same small candidates FASTA
        if not build_faidx(target_fasta, threads):
            sys.exit(1)

    # --- Step 3: check for duplicate query names ---
    if pairs_df["qname"].duplicated().any():
        dups = pairs_df[pairs_df["qname"].duplicated()]["qname"].unique()
        print(f"Error: duplicate qname values: {', '.join(map(str, dups))}", file=sys.stderr)
        sys.exit(1)

    # --- Step 4: parallel per-pair nucmer alignment ---
    print(f"Aligning {len(pairs_df)} pairs with {threads} workers...", file=sys.stderr)
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(threads, len(pairs_df))) as executor:
        futures = {
            executor.submit(
                pull_aln, idx, r["qname"], r["tname"],
                query_fasta, target_fasta, fastadir, alndir
            ): r["qname"]
            for idx, r in pairs_df.iterrows()
        }
        for future in concurrent.futures.as_completed(futures):
            qname = futures[future]
            try:
                rc = future.result()
                if rc != 0:
                    print(f"Alignment failed for {qname} (exit {rc})", file=sys.stderr)
            except Exception as e:
                print(f"Alignment raised exception for {qname}: {e}", file=sys.stderr)

    # --- Step 5: parse nucmer output → BED ---
    bed_df = parse_nucmer(alndir)
    if bed_df.empty:
        print("Error: no alignments produced. BED file will not be written.", file=sys.stderr)
        sys.exit(1)

    bed_path = os.path.join(outdir, "trimming.bed")
    bed_df.to_csv(bed_path, sep='\t', index=False, header=False)
    print(f"{len(bed_df)} BED intervals written to {bed_path}", file=sys.stderr)

    # --- Step 6: trim with seqkit subseq ---
    run_trimming(query_fasta, outdir, threads)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
