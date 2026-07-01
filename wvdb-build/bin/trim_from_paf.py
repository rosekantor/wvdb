#!/usr/bin/env python3
"""
trim_from_paf.py — select best trimming coordinates per cluster from a
minimap2 PAF and produce a trimmed FASTA via seqkit subseq.

For each cluster, considers all available pairwise alignments:
  rank1 vs rank2  (pairs12) → trim rank1
  rank1 vs rank3  (pairs13) → trim rank1
  rank2 vs rank3  (pairs23) → trim rank2

Picks the alignment with the longest alignment block length (PAF col 10).
Writes a BED file and runs seqkit subseq to produce trimmed.fasta.

Coordinate conversion:
  PAF: 0-based half-open [start, end)
  BED (seqkit subseq): 1-based closed [start, end]
  Conversion: bed_start = paf_start + 1,  bed_end = paf_end  (unchanged)

Usage:
  trim_from_paf.py \\
      --paf        alignments.paf \\
      --pairs12    pairs12.tsv \\
      --pairs13    pairs13.tsv    (optional) \\
      --pairs23    pairs23.tsv    (optional) \\
      --fasta      trimming_candidates.fasta \\
      --outdir     .
      --threads    8
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# PAF loading
# ---------------------------------------------------------------------------

PAF_COLS = [
    'qname', 'qlen', 'qstart', 'qend', 'strand',
    'tname', 'tlen', 'tstart', 'tend',
    'n_matches', 'aln_block_len', 'mapq'
]

def load_paf(paf_path):
    """Load minimap2 PAF into a DataFrame. Handles plain and bgzipped."""
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
    # Drop self-hits
    df = df[df['qname'] != df['tname']].reset_index(drop=True)
    return df


def load_pairs(tsv_path):
    """Return list of (query, target) tuples from a two-column TSV."""
    if not tsv_path:
        return []
    p = Path(tsv_path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    pairs = []
    with open(p) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                pairs.append((parts[0], parts[1]))
    return pairs


# ---------------------------------------------------------------------------
# BED building
# ---------------------------------------------------------------------------

def build_bed(pairs12, pairs13, pairs23, paf_df):
    """
    For each cluster, find the best PAF alignment across all pair types and
    return a BED DataFrame with columns [seq_id, bed_start, bed_end].

    Tiebreaker: longest alignment block length (PAF col 10).
    For rank12/rank13: trim the rank1 sequence (query coords used).
    For rank23:        trim the rank2 sequence (query coords used).

    PAF may report a pair in either orientation (query→target or
    target→query). We always extract query-side coordinates because
    the query is the sequence we intend to trim.
    """
    # Build lookup: (qname, tname) → best row by aln_block_len
    # Index both orientations so we find pairs regardless of PAF direction.
    best = {}   # (seq_to_trim, reference) → row
    for _, row in paf_df.iterrows():
        q, t = row['qname'], row['tname']
        for key, trim_seq, coords in [
            ((q, t), q, (row['qstart'], row['qend'])),
            ((t, q), t, (row['tstart'], row['tend'])),
        ]:
            if key not in best or row['aln_block_len'] > best[key][0]:
                best[key] = (row['aln_block_len'], trim_seq, coords[0], coords[1])

    # For each cluster, collect candidate alignments from all pair types
    # cluster_key → (aln_block_len, seq_id, paf_start, paf_end)
    candidates = {}   # cluster_rep (rank1) → best candidate so far

    def update(cluster_rep, seq_to_trim, ref):
        """Update candidates if this alignment is the best seen for the cluster."""
        key = (seq_to_trim, ref)
        if key not in best:
            return
        aln_len, trim_seq, paf_start, paf_end = best[key]
        if cluster_rep not in candidates or aln_len > candidates[cluster_rep][0]:
            candidates[cluster_rep] = (aln_len, trim_seq, paf_start, paf_end)

    # pairs12: rank1 (query) vs rank2 (target) — cluster_rep = rank1
    for rank1, rank2 in pairs12:
        update(rank1, rank1, rank2)

    # pairs13: rank1 (query) vs rank3 (target) — cluster_rep = rank1
    for rank1, rank3 in pairs13:
        update(rank1, rank1, rank3)

    # pairs23: rank2 (query) vs rank3 (target) — cluster_rep = rank1
    # need rank1→rank2 mapping to identify the cluster
    rank2_to_rank1 = {rank2: rank1 for rank1, rank2 in pairs12}
    for rank2, rank3 in pairs23:
        rank1 = rank2_to_rank1.get(rank2)
        if rank1:
            update(rank1, rank2, rank3)

    if not candidates:
        return pd.DataFrame(columns=['seq_id', 'bed_start', 'bed_end']), {}

    rows = []
    source_map = {}   # cluster_rep → seq_id actually trimmed
    n_rank1 = n_rank2 = 0
    for rank1, (aln_len, trim_seq, paf_start, paf_end) in candidates.items():
        bed_start = paf_start + 1   # 0-based → 1-based
        bed_end   = paf_end         # half-open → closed (end unchanged)
        rows.append({'seq_id': trim_seq, 'bed_start': bed_start, 'bed_end': bed_end})
        source_map[rank1] = trim_seq
        if trim_seq == rank1:
            n_rank1 += 1
        else:
            n_rank2 += 1

    print(f"trim_from_paf: {len(candidates)} clusters with alignments "
          f"({n_rank1} rank1 trimmed, {n_rank2} rank2 trimmed)", file=sys.stderr)

    no_hit = [r1 for r1 in {r1 for r1, _ in pairs12} if r1 not in candidates]
    if no_hit:
        print(f"trim_from_paf: {len(no_hit)} clusters had no PAF alignment "
              f"— these will be routed to blast_trim", file=sys.stderr)

    bed_df = pd.DataFrame(rows, columns=['seq_id', 'bed_start', 'bed_end'])
    return bed_df, source_map


# ---------------------------------------------------------------------------
# Trimming
# ---------------------------------------------------------------------------

def run_trimming(query_fasta, bed_path, outdir, threads):
    """Run seqkit subseq and clean coordinate suffixes from headers."""
    temp_fasta  = os.path.join(outdir, 'trimmed_temp.fasta')
    final_fasta = os.path.join(outdir, 'trimmed.fasta')

    cmd = ['seqkit', 'subseq', '--bed', bed_path, query_fasta,
           '-j', str(threads)]
    try:
        with open(temp_fasta, 'w') as fout:
            result = subprocess.run(cmd, stdout=fout,
                                    stderr=subprocess.PIPE, text=True)
    except OSError as e:
        print(f"Error running seqkit subseq: {e}", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(f"Error running seqkit subseq: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # seqkit subseq appends _start-end:. to headers — strip it
    pattern = re.compile(r'_\d+-\d+:\.')
    with open(temp_fasta) as fin, open(final_fasta, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                line = pattern.sub('', line)
            fout.write(line)
    os.remove(temp_fasta)
    return final_fasta


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Select best trimming per cluster from minimap2 PAF."
    )
    p.add_argument('--paf',     required=True, help='minimap2 PAF (plain or .gz)')
    p.add_argument('--pairs12', required=True, help='pairs12.tsv (rank1 vs rank2)')
    p.add_argument('--pairs13', default=None,  help='pairs13.tsv (rank1 vs rank3)')
    p.add_argument('--pairs23', default=None,  help='pairs23.tsv (rank2 vs rank3)')
    p.add_argument('--fasta',   required=True, help='trimming_candidates.fasta')
    p.add_argument('--outdir',  required=True, help='output directory')
    p.add_argument('--threads', type=int, default=1)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading PAF: {args.paf}", file=sys.stderr)
    paf_df = load_paf(args.paf)
    print(f"  {len(paf_df):,} non-self alignment records", file=sys.stderr)

    pairs12 = load_pairs(args.pairs12)
    pairs13 = load_pairs(args.pairs13)
    pairs23 = load_pairs(args.pairs23)
    print(f"  pairs12={len(pairs12)}, pairs13={len(pairs13)}, "
          f"pairs23={len(pairs23)}", file=sys.stderr)

    if not pairs12:
        print("No pairs to process — writing empty outputs.", file=sys.stderr)
        open(os.path.join(args.outdir, 'trimming.bed'),  'w').close()
        open(os.path.join(args.outdir, 'trimmed.fasta'), 'w').close()
        sys.exit(0)

    bed_df, source_map = build_bed(pairs12, pairs13, pairs23, paf_df)

    if bed_df.empty:
        print("Error: no PAF alignments found for any cluster pair. "
              "Check that the PAF was generated from the same FASTA.",
              file=sys.stderr)
        sys.exit(1)

    bed_path = os.path.join(args.outdir, 'trimming.bed')
    bed_df.to_csv(bed_path, sep='\t', index=False, header=False)
    print(f"  BED written: {len(bed_df)} intervals → {bed_path}", file=sys.stderr)

    run_trimming(args.fasta, bed_path, args.outdir, args.threads)
    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
