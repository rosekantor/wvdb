#!/usr/bin/env python3
"""
pick_best_trim1.py — select the best rank1 trimming per cluster.

For each cluster, compares trim_12 (rank1 trimmed against rank2) and trim_13
(rank1 trimmed against rank3, if available) and picks the longer result.
Then cross-references CheckV completeness to split into:

  complete_reps.fasta         → fed to reclustering (step 6)
  incomplete_rank1_ids.txt    → clusters needing trim_23 attempt
  branch_c_ids.txt            → 2-member clusters where trim_12 is incomplete
                                (no rank3 available; go directly to branch C)

Usage:
  pick_best_trim1.py \\
      --trim12-fasta    trim12/trimmed.fasta \\
      --trim13-fasta    trim13/trimmed.fasta \\
      --checkv-tsv      checkv_out/quality_summary.tsv \\
      --pairs12         pairs12.tsv \\
      --pairs13         pairs13.tsv \\
      --threshold       90 \\
      --complete-fasta  complete_reps.fasta \\
      --incomplete-ids  incomplete_rank1_ids.txt \\
      --branch-c-ids    branch_c_ids.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_fasta_lengths(fasta_path):
    """
    Return a dict {seq_id: length} for all records in fasta_path.
    Returns empty dict if fasta_path is None or does not exist.
    """
    if not fasta_path or not Path(fasta_path).exists():
        return {}
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, 'fasta')}


def load_fasta_records(fasta_path):
    """
    Return a dict {seq_id: SeqRecord} for all records in fasta_path.
    Returns empty dict if fasta_path is None or does not exist.
    """
    if not fasta_path or not Path(fasta_path).exists():
        return {}
    return {rec.id: rec for rec in SeqIO.parse(fasta_path, 'fasta')}


def load_pairs(tsv_path):
    """
    Load a two-column TSV (query, target) and return a dict {query: target}.
    Returns empty dict if tsv_path is None or does not exist.
    """
    if not tsv_path or not Path(tsv_path).exists():
        return {}
    pairs = {}
    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                pairs[parts[0]] = parts[1]   # query → target
    return pairs


def load_checkv_completeness(tsv_path, threshold):
    """
    Read CheckV quality_summary.tsv and return:
      complete_ids:   set of contig_ids with completeness > threshold
      incomplete_ids: set of contig_ids with completeness <= threshold or NA
    """
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)

    required = {'contig_id', 'completeness'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing columns in CheckV TSV: {sorted(missing)}. "
            f"Found: {list(df.columns)}"
        )

    comp = pd.to_numeric(df['completeness'], errors='coerce')
    ids  = df['contig_id'].astype(str)

    complete_ids   = set(ids[comp.notna() & (comp > threshold)])
    incomplete_ids = set(ids[comp.isna()  | (comp <= threshold)])

    return complete_ids, incomplete_ids


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def pick_best(rank1_id, trim12_lengths, trim13_lengths):
    """
    For a given rank1_id, return the ID of whichever trimming is longer.
    If only trim12 is available, return rank1_id (the trim12 result carries
    the same ID as the original rank1 sequence after name cleanup).
    Returns (best_id, source) where source is 'trim12' or 'trim13'.
    """
    len12 = trim12_lengths.get(rank1_id)
    len13 = trim13_lengths.get(rank1_id)

    if len12 is None and len13 is None:
        return None, None
    if len12 is None:
        return rank1_id, 'trim13'
    if len13 is None:
        return rank1_id, 'trim12'

    # Both available: prefer longer; tie goes to trim12
    if len13 > len12:
        return rank1_id, 'trim13'
    return rank1_id, 'trim12'


def main():
    parser = argparse.ArgumentParser(
        description="Select best rank1 trimming per cluster and split by CheckV completeness."
    )
    parser.add_argument('--trim12-fasta',   required=True,
                        help='trimmed.fasta from TRIM_GENOMES_12')
    parser.add_argument('--trim13-fasta',   default=None,
                        help='trimmed.fasta from TRIM_GENOMES_13 (omit if no 3+-member clusters)')
    parser.add_argument('--checkv-tsv',     required=True,
                        help='CheckV quality_summary.tsv for the merged trim1 candidates')
    parser.add_argument('--pairs12',        required=True,
                        help='pairs12.tsv (query=rank1, target=rank2)')
    parser.add_argument('--pairs13',        default=None,
                        help='pairs13.tsv (query=rank1, target=rank3)')
    parser.add_argument('--threshold',      type=float, default=90.0,
                        help='CheckV completeness threshold (default: 90)')
    parser.add_argument('--complete-fasta', required=True,
                        help='Output FASTA of complete rank1 trimmed representatives')
    parser.add_argument('--incomplete-ids', required=True,
                        help='Output file: rank1 IDs whose best trimming is incomplete '
                             'AND have a rank3 available (→ try trim23)')
    parser.add_argument('--branch-c-ids',   required=True,
                        help='Output file: rank1 IDs whose trim12 is incomplete '
                             'AND no rank3 available (→ branch C directly)')
    args = parser.parse_args()

    # --- Load inputs ---
    pairs12      = load_pairs(args.pairs12)
    pairs13      = load_pairs(args.pairs13)
    trim12_lens  = load_fasta_lengths(args.trim12_fasta)
    trim13_lens  = load_fasta_lengths(args.trim13_fasta)
    trim12_recs  = load_fasta_records(args.trim12_fasta)
    trim13_recs  = load_fasta_records(args.trim13_fasta)

    complete_ids, incomplete_ids = load_checkv_completeness(args.checkv_tsv, args.threshold)

    # rank1 IDs that have a rank3 (i.e. appear in pairs13)
    has_rank3 = set(pairs13.keys())

    # All rank1 IDs seen across both pair files
    all_rank1_ids = set(pairs12.keys()) | set(pairs13.keys())

    # --- Per-cluster decision ---
    complete_records   = []
    incomplete_need23  = []   # → incomplete_rank1_ids.txt
    incomplete_branch_c = []  # → branch_c_ids.txt

    n_trim12_wins = 0
    n_trim13_wins = 0
    n_no_aln      = 0

    for rank1_id in sorted(all_rank1_ids):
        best_id, source = pick_best(rank1_id, trim12_lens, trim13_lens)

        if best_id is None:
            # No trimming succeeded for this rank1 at all
            n_no_aln += 1
            print(
                f"Warning: no trimmed sequence found for {rank1_id} — "
                f"sending to branch C.",
                file=sys.stderr
            )
            incomplete_branch_c.append(rank1_id)
            continue

        if source == 'trim12':
            n_trim12_wins += 1
            rec = trim12_recs.get(best_id)
        else:
            n_trim13_wins += 1
            rec = trim13_recs.get(best_id)

        if rec is None:
            print(
                f"Warning: {rank1_id} selected from {source} but record missing — "
                f"sending to branch C.",
                file=sys.stderr
            )
            incomplete_branch_c.append(rank1_id)
            continue

        if best_id in complete_ids:
            complete_records.append(rec)
        else:
            # Incomplete — decide where to send based on rank3 availability
            if rank1_id in has_rank3:
                incomplete_need23.append(rank1_id)
            else:
                incomplete_branch_c.append(rank1_id)

    # --- Write outputs ---
    Path(args.complete_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(args.complete_fasta, 'w') as f:
        SeqIO.write(complete_records, f, 'fasta')

    with open(args.incomplete_ids, 'w') as f:
        for id_ in sorted(incomplete_need23):
            f.write(id_ + '\n')

    with open(args.branch_c_ids, 'w') as f:
        for id_ in sorted(incomplete_branch_c):
            f.write(id_ + '\n')

    # --- Summary ---
    print(f"[pick_best_trim1] threshold: {args.threshold}", file=sys.stderr)
    print(f"[pick_best_trim1] rank1 IDs processed: {len(all_rank1_ids)}", file=sys.stderr)
    print(f"[pick_best_trim1] trim12 selected: {n_trim12_wins}", file=sys.stderr)
    print(f"[pick_best_trim1] trim13 selected: {n_trim13_wins}", file=sys.stderr)
    print(f"[pick_best_trim1] no alignment found: {n_no_aln}", file=sys.stderr)
    print(f"[pick_best_trim1] complete reps written: {len(complete_records)}", file=sys.stderr)
    print(f"[pick_best_trim1] incomplete → try trim23: {len(incomplete_need23)}", file=sys.stderr)
    print(f"[pick_best_trim1] incomplete → branch C: {len(incomplete_branch_c)}", file=sys.stderr)

    if len(complete_ids) > 0 and len(complete_records) == 0:
        raise RuntimeError(
            "CheckV reported complete IDs but zero records were written. "
            "Possible ID mismatch between trimmed FASTA headers and CheckV contig_id column."
        )


if __name__ == '__main__':
    main()
