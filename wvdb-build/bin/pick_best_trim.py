#!/usr/bin/env python3
"""
pick_best_trim.py — split cluster trim output by CheckV completeness.

All pairs (rank12, rank13, rank23) were already evaluated by trim_from_paf.py
which chose the best alignment per cluster. This script simply reads the
CheckV result and routes each cluster to:

  complete_reps.fasta   → step 5 (recluster)
  blast_trim_ids.txt    → step 4 (blast_trim): incomplete or no alignment

Usage:
  pick_best_trim.py \\
      --trimmed-fasta  trimmed.fasta \\
      --checkv-tsv     checkv_out/quality_summary.tsv \\
      --threshold      90 \\
      --complete-fasta complete_reps.fasta \\
      --blast-trim-ids blast_trim_ids.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--trimmed-fasta',  required=True)
    p.add_argument('--checkv-tsv',     required=True)
    p.add_argument('--threshold',      type=float, default=90.0)
    p.add_argument('--complete-fasta', required=True)
    p.add_argument('--blast-trim-ids', required=True)
    args = p.parse_args()

    # Load trimmed sequences
    trimmed = {r.id: r for r in SeqIO.parse(args.trimmed_fasta, 'fasta')}

    # Load CheckV completeness
    df = pd.read_csv(args.checkv_tsv, sep='\t', dtype=str)
    comp = pd.to_numeric(df['completeness'], errors='coerce')
    ids  = df['contig_id'].astype(str)
    complete_ids   = set(ids[comp.notna() & (comp > args.threshold)])
    incomplete_ids = set(ids[comp.isna()  | (comp <= args.threshold)])

    complete_recs = []
    blast_trim    = []

    for seq_id, rec in trimmed.items():
        if seq_id in complete_ids:
            complete_recs.append(rec)
        else:
            blast_trim.append(seq_id)

    with open(args.complete_fasta, 'w') as f:
        SeqIO.write(complete_recs, f, 'fasta')
    with open(args.blast_trim_ids, 'w') as f:
        f.writelines(id_ + '\n' for id_ in sorted(blast_trim))

    print(f"[pick_best_trim] threshold       : {args.threshold}", file=sys.stderr)
    print(f"[pick_best_trim] trimmed seqs    : {len(trimmed)}", file=sys.stderr)
    print(f"[pick_best_trim] complete        : {len(complete_recs)}", file=sys.stderr)
    print(f"[pick_best_trim] → blast_trim    : {len(blast_trim)}", file=sys.stderr)

    if complete_ids and not complete_recs:
        raise RuntimeError(
            "CheckV reported complete IDs but zero records written. "
            "Possible ID mismatch between trimmed FASTA and CheckV contig_id."
        )


if __name__ == '__main__':
    main()
