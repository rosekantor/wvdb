#!/usr/bin/env python3
"""
pick_best_trim.py — select the best rank1 trimming per cluster from
trim12, trim13, and trim23 results, using CheckV completeness as the
primary criterion and sequence length as a tiebreaker.

The combined FASTA fed to CheckV uses suffixed IDs to distinguish the
three trimmings of the same rank1 sequence:
  <rank1_id>_trim12
  <rank1_id>_trim13
  <rank1_id>_trim23

These suffixes are stripped from the winning sequence before writing output.

Selection priority per rank1 cluster:
  1. Higher CheckV completeness wins
  2. Equal completeness → longer sequence wins
  3. Complete tie → prefer trim12 > trim13 > trim23

Outputs:
  complete_reps.fasta    sequences above completeness threshold → step 5
  blast_trim_ids.txt     rank1 IDs below threshold → step 4 (blast trim)

Usage:
  pick_best_trim.py \\
      --trim12    trim12.trimmed.fasta \\
      --trim13    trim13.trimmed.fasta  (optional) \\
      --trim23    trim23.trimmed.fasta  (optional) \\
      --checkv    checkv_out/quality_summary.tsv \\
      --threshold 90 \\
      --complete  complete_reps.fasta \\
      --blast-ids blast_trim_ids.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def load_fasta(path):
    """Return {id: SeqRecord}. Empty dict if path is None/missing/empty."""
    if not path:
        return {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    return {rec.id: rec for rec in SeqIO.parse(path, 'fasta')}


def load_pairs(path):
    """Return {rank1: rank2} from a two-column TSV. Empty dict if missing."""
    if not path:
        return {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    pairs = {}
    with open(p) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                pairs[parts[0]] = parts[1]
    return pairs


def load_checkv(path):
    """Return {contig_id: completeness_float}. NA/missing → -1.0."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    comp = pd.to_numeric(df['completeness'], errors='coerce').fillna(-1.0)
    return dict(zip(df['contig_id'].astype(str), comp))
    """Return {contig_id: completeness_float}. NA/missing → -1.0."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    comp = pd.to_numeric(df['completeness'], errors='coerce').fillna(-1.0)
    return dict(zip(df['contig_id'].astype(str), comp))


def best_of(candidates):
    """
    Select the best (rec, source, comp) from a list of (rec, source, comp)
    tuples where rec is not None.
    Priority: highest comp → longest seq → earliest in list (trim12>13>23).
    """
    best = None
    for rec, source, comp in candidates:
        if rec is None:
            continue
        if best is None:
            best = (rec, source, comp)
            continue
        b_rec, b_source, b_comp = best
        if comp > b_comp:
            best = (rec, source, comp)
        elif comp == b_comp and len(rec.seq) > len(b_rec.seq):
            best = (rec, source, comp)
    return best


def main():
    p = argparse.ArgumentParser(
        description="Select best rank1 trimming from trim12/13/23 by CheckV completeness."
    )
    p.add_argument('--trim12',    required=True,  help='trim12.trimmed.fasta')
    p.add_argument('--trim13',    default=None,   help='trim13.trimmed.fasta (optional)')
    p.add_argument('--trim23',    default=None,   help='trim23.trimmed.fasta (optional)')
    p.add_argument('--checkv',    required=True,  help='CheckV quality_summary.tsv')
    p.add_argument('--pairs12',   default=None,   help='pairs12.tsv for trim23 rank2→rank1 lookup')
    p.add_argument('--trim-log',  default=None,   help='Output: TSV log of trim source per rank1 ID')
    p.add_argument('--threshold', type=float, default=90.0,
                   help='CheckV completeness threshold (default: 90)')
    p.add_argument('--complete',  required=True,  help='Output: complete_reps.fasta')
    p.add_argument('--blast-ids', required=True,  help='Output: blast_trim_ids.txt')
    args = p.parse_args()

    trim12 = load_fasta(args.trim12)
    trim13 = load_fasta(args.trim13)
    trim23 = load_fasta(args.trim23)  # keyed by rank2 IDs
    checkv = load_checkv(args.checkv)

    # Build rank1 → rank2 map so we can look up trim23 sequences by rank2 ID
    pairs12       = load_pairs(args.pairs12)
    rank1_to_rank2 = pairs12  # {rank1: rank2}

    # all_rank1 comes from trim12 and trim13 only — trim23 uses rank2 IDs
    # which have already been translated to rank1 in merge_trimmed.py via
    # the _trim23 suffix on rank1 IDs in all_trimmed.fasta → CheckV output.
    # The trim23 FASTA here still has rank2 IDs; we only use it for sequence
    # lookup after the rank1 ID is established from trim12/trim13.
    all_rank1 = set(trim12) | set(trim13)

    wins      = {'trim12': 0, 'trim13': 0, 'trim23': 0}
    n_complete   = 0
    n_blast_trim = 0
    n_no_aln     = 0

    complete_records = []
    blast_trim_ids   = []
    log_rows         = []

    for rank1_id in sorted(all_rank1):
        rec12 = trim12.get(rank1_id)
        rec13 = trim13.get(rank1_id)
        # trim23 is keyed by rank2 ID — look up via pairs12 map
        rank2_id = rank1_to_rank2.get(rank1_id)
        rec23 = trim23.get(rank2_id) if rank2_id else None

        if rec12 is None and rec13 is None and rec23 is None:
            n_no_aln += 1
            print(f"Warning: no trimmed sequence for {rank1_id} — sending to blast_trim",
                  file=sys.stderr)
            blast_trim_ids.append(rank1_id)
            continue

        # Look up CheckV completeness for each available trimming
        comp12 = checkv.get(f"{rank1_id}_trim12", -1.0) if rec12 else -1.0
        comp13 = checkv.get(f"{rank1_id}_trim13", -1.0) if rec13 else -1.0
        comp23 = checkv.get(f"{rank1_id}_trim23", -1.0) if rec23 else -1.0

        # Select best — order determines tie preference (trim12 > trim13 > trim23)
        result = best_of([
            (rec12, 'trim12', comp12),
            (rec13, 'trim13', comp13),
            (rec23, 'trim23', comp23),
        ])

        best_rec, source, best_comp = result
        wins[source] += 1

        # Restore original (non-suffixed) ID
        best_rec = best_rec.__class__(best_rec.seq, id=rank1_id, description='')

        if best_comp >= args.threshold:
            complete_records.append(best_rec)
            n_complete += 1
            outcome = 'complete'
        else:
            blast_trim_ids.append(rank1_id)
            n_blast_trim += 1
            outcome = 'blast_trim'

        log_rows.append({
            'rank1_id':    rank1_id,
            'source':      source,
            'completeness': round(best_comp, 2),
            'seq_len':     len(best_rec.seq),
            'outcome':     outcome,
        })

    with open(args.complete, 'w') as f:
        SeqIO.write(complete_records, f, 'fasta')
    with open(args.blast_ids, 'w') as f:
        f.writelines(id_ + '\n' for id_ in sorted(blast_trim_ids))

    if args.trim_log:
        import csv
        with open(args.trim_log, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['rank1_id', 'source',
                                                    'completeness', 'seq_len', 'outcome'],
                                    delimiter='\t')
            writer.writeheader()
            writer.writerows(log_rows)

    print(f"[pick_best_trim] rank1 IDs processed : {len(all_rank1)}", file=sys.stderr)
    print(f"[pick_best_trim] trim12 selected      : {wins['trim12']}", file=sys.stderr)
    print(f"[pick_best_trim] trim13 selected      : {wins['trim13']}", file=sys.stderr)
    print(f"[pick_best_trim] trim23 selected      : {wins['trim23']}", file=sys.stderr)
    print(f"[pick_best_trim] no alignment         : {n_no_aln}", file=sys.stderr)
    print(f"[pick_best_trim] complete → recluster : {n_complete}", file=sys.stderr)
    print(f"[pick_best_trim] incomplete → blast   : {n_blast_trim}", file=sys.stderr)


if __name__ == '__main__':
    main()
