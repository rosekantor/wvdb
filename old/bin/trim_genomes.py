#!/usr/bin/env python3
"""
trim_genomes.py — trim query sequences using pairwise alignments from a minimap2 PAF file.

Replaces the nucmer/show-coords based approach with PAF parsing, eliminating per-pair
sequence extraction and enabling parallel trimming from a single shared alignment file.

Usage:
  trim_genomes.py -p alignments.paf -pairs pairs12.tsv -f candidates.fasta -o outdir [-t THREADS]

Input PAF format (minimap2 -x asm5 output):
  Col 0:  query name
  Col 1:  query length
  Col 2:  query start (0-based)
  Col 3:  query end   (0-based, open)
  Col 4:  strand (+/-)
  Col 5:  target name
  Col 6:  target length
  Col 7:  target start (0-based)
  Col 8:  target end   (0-based, open)
  Col 9:  residue matches
  Col 10: alignment block length
  Col 11: mapping quality

Coordinate conversion:
  PAF uses 0-based half-open intervals [start, end).
  seqkit subseq --bed uses 1-based closed intervals [start, end].
  Conversion: bed_start = paf_start + 1,  bed_end = paf_end  (no change to end)

Trimming logic:
  For each (query, target) pair in pairs_tsv:
    - Find all PAF rows where col0==query AND col5==target (or reverse orientation)
    - Keep the row with the largest alignment block length (col10)
    - The BED interval is derived from the QUERY coordinates (we are trimming the query)
    - Write BED: query_name <TAB> bed_start <TAB> bed_end
  Run seqkit subseq --bed on query_fasta to produce trimmed.fasta.
  Clean up coordinate suffixes from FASTA headers added by seqkit subseq.
"""

import argparse
import os
import re
import subprocess
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# PAF parsing
# ---------------------------------------------------------------------------

PAF_COLS = [
    'qname', 'qlen', 'qstart', 'qend', 'strand',
    'tname', 'tlen', 'tstart', 'tend',
    'n_matches', 'aln_block_len', 'mapq'
]

def load_paf(paf_path):
    """
    Load a minimap2 PAF file into a DataFrame.
    Only the first 12 columns are read; optional cs/cg tags are ignored.
    Handles bgzipped PAF (.paf.gz) transparently via pandas.
    """
    df = pd.read_csv(
        paf_path,
        sep='\t',
        header=None,
        usecols=range(12),
        names=PAF_COLS,
        dtype={
            'qname': str, 'tname': str, 'strand': str,
            'qlen': int, 'qstart': int, 'qend': int,
            'tlen': int, 'tstart': int, 'tend': int,
            'n_matches': int, 'aln_block_len': int, 'mapq': int,
        }
    )
    return df


def load_pairs(pairs_tsv):
    """
    Load a two-column TSV of (query, target) pairs.
    Returns a list of (query, target) tuples.
    """
    pairs = []
    with open(pairs_tsv) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                print(f"Warning: skipping malformed pairs line: {repr(line)}", file=sys.stderr)
                continue
            pairs.append((parts[0], parts[1]))
    return pairs


def build_bed(pairs, paf_df):
    """
    For each (query, target) pair, find the best PAF alignment and derive
    a BED interval on the query sequence.

    PAF is directional: minimap2 reports the alignment from query→target.
    We always want to trim the query, so we use query coordinates regardless
    of strand. This is correct because PAF query coordinates are always on
    the forward strand of the query sequence.

    Returns a DataFrame with columns [qname, bed_start, bed_end] and a list
    of (query, target) pairs that had no PAF hit.
    """
    # Build a lookup dict: (qname, tname) → best row (max aln_block_len)
    # Also index (tname, qname) in case minimap2 reported the pair reversed.
    paf_indexed = {}
    for _, row in paf_df.iterrows():
        key_fwd = (row['qname'], row['tname'])
        key_rev = (row['tname'], row['qname'])
        for key in (key_fwd, key_rev):
            if key not in paf_indexed or row['aln_block_len'] > paf_indexed[key]['aln_block_len']:
                paf_indexed[key] = row

    bed_rows = []
    no_hit = []

    for query, target in pairs:
        key = (query, target)
        row = paf_indexed.get(key)

        if row is None:
            no_hit.append((query, target))
            continue

        # Determine which sequence is the query in this PAF row
        if row['qname'] == query:
            paf_start = row['qstart']
            paf_end   = row['qend']
        else:
            # Pair was found reversed in PAF; use target coords as query coords
            # (minimap2 reported target→query; we still trim 'query')
            paf_start = row['tstart']
            paf_end   = row['tend']

        # Convert 0-based half-open [paf_start, paf_end) → 1-based closed [bed_start, bed_end]
        bed_start = paf_start + 1
        bed_end   = paf_end       # end is the same in 1-based closed

        bed_rows.append({'qname': query, 'bed_start': bed_start, 'bed_end': bed_end})

    if no_hit:
        print(
            f"Warning: {len(no_hit)} pairs had no PAF alignment and will be skipped:\n" +
            '\n'.join(f"  {q}\t{t}" for q, t in no_hit[:10]) +
            ('\n  ...' if len(no_hit) > 10 else ''),
            file=sys.stderr
        )

    bed_df = pd.DataFrame(bed_rows, columns=['qname', 'bed_start', 'bed_end'])
    return bed_df, no_hit


# ---------------------------------------------------------------------------
# Trimming
# ---------------------------------------------------------------------------

def run_trimming(query_fasta, bed_path, outdir, threads):
    """
    Run seqkit subseq with the BED file to trim sequences, then clean up
    coordinate suffixes added by seqkit to the FASTA headers.
    """
    temp_fasta  = os.path.join(outdir, 'trimmed_temp.fasta')
    final_fasta = os.path.join(outdir, 'trimmed.fasta')

    cmd_subseq = [
        'seqkit', 'subseq',
        '--bed', bed_path,
        query_fasta,
        '-j', str(threads),
    ]

    try:
        with open(temp_fasta, 'w') as fout:
            result = subprocess.run(
                cmd_subseq,
                stdout=fout,
                stderr=subprocess.PIPE,
                text=True,
            )
    except OSError as e:
        print(f"Error writing temporary trimmed FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(f"Error running seqkit subseq: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # seqkit subseq appends '_start-end:.' to each header; strip it
    pattern = re.compile(r'_\d+-\d+:\.')

    try:
        with open(temp_fasta) as fin, open(final_fasta, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    line = pattern.sub('', line)
                fout.write(line)
    except OSError as e:
        print(f"Error rewriting trimmed FASTA: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        os.remove(temp_fasta)
    except OSError:
        print(f"Warning: could not remove temp file '{temp_fasta}'", file=sys.stderr)

    return final_fasta


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Trim query sequences using pairwise alignments from a minimap2 PAF file. "
            "Replaces the nucmer/show-coords approach used in the original pipeline."
        )
    )
    parser.add_argument(
        '-p', '--paf', required=True,
        help='minimap2 PAF file (plain or bgzipped .paf.gz)'
    )
    parser.add_argument(
        '--pairs', required=True,
        help='Two-column TSV of (query, target) pairs to trim'
    )
    parser.add_argument(
        '-f', '--fasta', required=True,
        help='FASTA containing the query sequences (candidates.fasta from seqkit grep)'
    )
    parser.add_argument(
        '-o', '--outdir', required=True,
        help='Output directory; trimmed.fasta and trimming.bed are written here'
    )
    parser.add_argument(
        '-t', '--threads', type=int, default=1,
        help='Threads to pass to seqkit subseq (default: 1)'
    )
    return parser.parse_args()


def main():
    args = parse_args()

    for path in [args.paf, args.pairs, args.fasta]:
        if not os.path.isfile(path):
            print(f"Error: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load PAF and pairs
    print(f"Loading PAF: {args.paf}", file=sys.stderr)
    paf_df = load_paf(args.paf)
    print(f"  {len(paf_df):,} alignment records", file=sys.stderr)

    pairs = load_pairs(args.pairs)
    print(f"  {len(pairs):,} pairs to trim", file=sys.stderr)

    if not pairs:
        print("No pairs to process — writing empty outputs.", file=sys.stderr)
        open(os.path.join(args.outdir, 'trimming.bed'),  'w').close()
        open(os.path.join(args.outdir, 'trimmed.fasta'), 'w').close()
        sys.exit(0)

    # 2. Build BED from PAF
    bed_df, no_hit = build_bed(pairs, paf_df)

    if bed_df.empty:
        print(
            "Error: no PAF alignments found for any pair. "
            "Check that the PAF was generated from the same FASTA as the pairs TSV.",
            file=sys.stderr
        )
        sys.exit(1)

    bed_path = os.path.join(args.outdir, 'trimming.bed')
    bed_df.to_csv(bed_path, sep='\t', index=False, header=False)
    print(f"  BED written: {len(bed_df)} intervals → {bed_path}", file=sys.stderr)

    # 3. Trim with seqkit subseq
    final_fasta = run_trimming(args.fasta, bed_path, args.outdir, args.threads)
    print(f"  Trimmed FASTA written → {final_fasta}", file=sys.stderr)


if __name__ == '__main__':
    main()
