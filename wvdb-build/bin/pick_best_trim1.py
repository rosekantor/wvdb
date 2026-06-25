#!/usr/bin/env python3
"""
pick_best_trim1.py — select the best rank1 trimming per cluster.

Compares trim12 (rank1 vs rank2) and trim13 (rank1 vs rank3) using CheckV
completeness as the primary criterion, then sequence length as a tiebreaker.

The combined FASTA fed to CheckV uses suffixed IDs to distinguish the two
trimmings of the same rank1 sequence:
  <rank1_id>_trim12
  <rank1_id>_trim13

Selection logic per rank1 cluster:
  1. Higher CheckV completeness wins
  2. Tie → longer sequence wins
  3. Tie → trim12 preferred
  If only trim12 available (2-member cluster): use trim12 directly

Outputs:
  complete_reps.fasta         → step 5 (recluster)
  incomplete_rank1_ids.txt    → FILTER_PAIRS23 → branch B (try trim23)
  blast_trim_ids.txt          → FETCH_BLAST_INPUT_SEQS → step 4 (blast trim)

Usage:
  pick_best_trim1.py \\
      --trim12-fasta    trim12.trimmed.fasta \\
      --trim13-fasta    trim13.trimmed.fasta \\
      --checkv-tsv      checkv_out/quality_summary.tsv \\
      --pairs12         pairs12.tsv \\
      --pairs13         pairs13.tsv \\
      --threshold       90 \\
      --complete-fasta  complete_reps.fasta \\
      --incomplete-ids  incomplete_rank1_ids.txt \\
      --blast-trim-ids  blast_trim_ids.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_fasta(fasta_path):
    """Return {id: SeqRecord} for all records. Empty dict if missing/empty."""
    if not fasta_path:
        return {}
    p = Path(fasta_path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    return {rec.id: rec for rec in SeqIO.parse(fasta_path, 'fasta')}


def load_pairs(tsv_path):
    """Return {query: target} from a two-column TSV. Empty dict if missing."""
    if not tsv_path:
        return {}
    p = Path(tsv_path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    pairs = {}
    with open(p) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                pairs[parts[0]] = parts[1]
    return pairs


def load_checkv(tsv_path):
    """
    Return {contig_id: completeness_float} from CheckV quality_summary.tsv.
    NA completeness stored as -1.0 so numeric comparisons always work.
    """
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)
    required = {'contig_id', 'completeness'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in CheckV TSV: {sorted(missing)}")
    comp = pd.to_numeric(df['completeness'], errors='coerce').fillna(-1.0)
    return dict(zip(df['contig_id'].astype(str), comp))


def build_combined_fasta(trim12_recs, trim13_recs, out_path):
    """
    Write a combined FASTA with suffixed IDs:
      <rank1_id>_trim12  (from trim12)
      <rank1_id>_trim13  (from trim13)
    This is the FASTA passed to CheckV so both trimmings are assessed together.
    Returns the set of rank1 IDs that have a trim13 record.
    """
    records = []
    has_trim13 = set()
    for rank1_id, rec in trim12_recs.items():
        r = rec.__class__(rec.seq, id=f"{rank1_id}_trim12",
                          description='')
        records.append(r)
    for rank1_id, rec in trim13_recs.items():
        r = rec.__class__(rec.seq, id=f"{rank1_id}_trim13",
                          description='')
        records.append(r)
        has_trim13.add(rank1_id)
    with open(out_path, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    return has_trim13


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Select best rank1 trimming per cluster by CheckV completeness."
    )
    parser.add_argument('--trim12-fasta',   required=True)
    parser.add_argument('--trim13-fasta',   default=None)
    parser.add_argument('--checkv-tsv',     required=True)
    parser.add_argument('--pairs12',        required=True)
    parser.add_argument('--pairs13',        default=None)
    parser.add_argument('--threshold',      type=float, default=90.0)
    parser.add_argument('--complete-fasta', required=True)
    parser.add_argument('--incomplete-ids', required=True)
    parser.add_argument('--blast-trim-ids', required=True)
    args = parser.parse_args()

    pairs12     = load_pairs(args.pairs12)
    pairs13     = load_pairs(args.pairs13)
    trim12_recs = load_fasta(args.trim12_fasta)
    trim13_recs = load_fasta(args.trim13_fasta)
    checkv      = load_checkv(args.checkv_tsv)

    has_rank3   = set(pairs13.keys())
    all_rank1   = set(pairs12.keys()) | has_rank3

    n_trim12_wins = 0
    n_trim13_wins = 0
    n_complete    = 0
    n_need23      = 0
    n_blast_trim  = 0
    n_no_aln      = 0

    complete_records  = []
    incomplete_ids    = []
    blast_trim_ids    = []

    for rank1_id in sorted(all_rank1):
        rec12 = trim12_recs.get(rank1_id)
        rec13 = trim13_recs.get(rank1_id)

        if rec12 is None and rec13 is None:
            n_no_aln += 1
            print(f"Warning: no trimmed sequence for {rank1_id} — sending to blast_trim",
                  file=sys.stderr)
            blast_trim_ids.append(rank1_id)
            continue

        # Look up CheckV completeness using suffixed IDs
        comp12 = checkv.get(f"{rank1_id}_trim12", -1.0) if rec12 else -1.0
        comp13 = checkv.get(f"{rank1_id}_trim13", -1.0) if rec13 else -1.0

        # Select best: completeness first, then length, then prefer trim12
        if rec12 is None:
            best_rec, source, best_comp = rec13, 'trim13', comp13
        elif rec13 is None:
            best_rec, source, best_comp = rec12, 'trim12', comp12
        elif comp13 > comp12:
            best_rec, source, best_comp = rec13, 'trim13', comp13
        elif comp12 > comp13:
            best_rec, source, best_comp = rec12, 'trim12', comp12
        elif len(rec13.seq) > len(rec12.seq):
            best_rec, source, best_comp = rec13, 'trim13', comp13
        else:
            best_rec, source, best_comp = rec12, 'trim12', comp12

        if source == 'trim12':
            n_trim12_wins += 1
        else:
            n_trim13_wins += 1

        # Restore original (non-suffixed) ID for output
        best_rec = best_rec.__class__(best_rec.seq, id=rank1_id, description='')

        if best_comp > args.threshold:
            complete_records.append(best_rec)
            n_complete += 1
        else:
            if rank1_id in has_rank3:
                incomplete_ids.append(rank1_id)
                n_need23 += 1
            else:
                blast_trim_ids.append(rank1_id)
                n_blast_trim += 1

    # Write outputs
    with open(args.complete_fasta, 'w') as f:
        SeqIO.write(complete_records, f, 'fasta')
    with open(args.incomplete_ids, 'w') as f:
        f.writelines(id_ + '\n' for id_ in sorted(incomplete_ids))
    with open(args.blast_trim_ids, 'w') as f:
        f.writelines(id_ + '\n' for id_ in sorted(blast_trim_ids))

    print(f"[pick_best_trim1] rank1 IDs processed : {len(all_rank1)}", file=sys.stderr)
    print(f"[pick_best_trim1] trim12 selected      : {n_trim12_wins}", file=sys.stderr)
    print(f"[pick_best_trim1] trim13 selected      : {n_trim13_wins}", file=sys.stderr)
    print(f"[pick_best_trim1] no alignment         : {n_no_aln}", file=sys.stderr)
    print(f"[pick_best_trim1] complete → recluster : {n_complete}", file=sys.stderr)
    print(f"[pick_best_trim1] incomplete → trim23  : {n_need23}", file=sys.stderr)
    print(f"[pick_best_trim1] incomplete → blast   : {n_blast_trim}", file=sys.stderr)

    if complete_records == [] and n_no_aln < len(all_rank1):
        print(
            "Warning: zero complete reps written but some trimmings existed. "
            "Check CheckV completeness values and --threshold.",
            file=sys.stderr
        )


if __name__ == '__main__':
    main()
