#!/usr/bin/env python3
"""
merge_trimmed.py — combine trim12, trim13, trim23 FASTAs into a single
FASTA with suffixed IDs for input to CheckV.

trim12 and trim13 contain rank1 IDs directly.
trim23 contains rank2 IDs — these are translated back to rank1 IDs using
pairs12.tsv (which maps rank1 → rank2) before suffixing.

Output IDs:
  <rank1_id>_trim12
  <rank1_id>_trim13
  <rank1_id>_trim23

This allows CheckV completeness scores to be attributed back to the
specific trimming that produced them in pick_best_trim.py, all keyed
on rank1 IDs.

Usage:
  merge_trimmed.py \\
      --trim12   trim12.trimmed.fasta \\
      --trim13   trim13.trimmed.fasta   (optional) \\
      --trim23   trim23.trimmed.fasta   (optional) \\
      --pairs12  pairs12.tsv            (required if --trim23 provided) \\
      --out      all_trimmed.fasta
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO


def load_pairs(path):
    """Return {rank1: rank2} from pairs TSV."""
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


def load_and_suffix(fasta_path, suffix, id_map=None):
    """
    Load FASTA, append suffix to each ID, return list of SeqRecords.
    If id_map is provided, translate IDs through it first
    (used to map rank2 → rank1 for trim23).
    Records whose ID is not in id_map are skipped with a warning.
    """
    if not fasta_path:
        return []
    p = Path(fasta_path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    recs = []
    skipped = 0
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        if id_map is not None:
            new_id = id_map.get(rec.id)
            if new_id is None:
                skipped += 1
                continue
            rec.id = new_id
        rec.id          = f"{rec.id}_{suffix}"
        rec.description = ''
        recs.append(rec)
    if skipped:
        print(f"[merge_trimmed] {suffix}: {skipped} records skipped (ID not in map)",
              file=sys.stderr)
    return recs


def main():
    p = argparse.ArgumentParser(
        description="Combine trim12/13/23 FASTAs with rank1-suffixed IDs for CheckV."
    )
    p.add_argument('--trim12',  required=True)
    p.add_argument('--trim13',  default=None)
    p.add_argument('--trim23',  default=None)
    p.add_argument('--pairs12', default=None,
                   help='pairs12.tsv (rank1→rank2); required to translate trim23 IDs')
    p.add_argument('--out',     required=True)
    args = p.parse_args()

    # Build rank2 → rank1 lookup from pairs12 (inverse of rank1 → rank2)
    rank2_to_rank1 = None
    if args.trim23 and not args.pairs12:
        print("Warning: --trim23 provided without --pairs12; trim23 IDs cannot be "
              "translated to rank1 and will be skipped", file=sys.stderr)
    elif args.pairs12:
        pairs12 = load_pairs(args.pairs12)
        rank2_to_rank1 = {v: k for k, v in pairs12.items()}

    recs  = load_and_suffix(args.trim12, 'trim12')
    recs += load_and_suffix(args.trim13, 'trim13')
    recs += load_and_suffix(args.trim23, 'trim23', id_map=rank2_to_rank1)

    if not recs:
        print("ERROR: no trimmed sequences found in any input FASTA", file=sys.stderr)
        sys.exit(1)

    with open(args.out, 'w') as f:
        SeqIO.write(recs, f, 'fasta')

    print(f"[merge_trimmed] {len(recs)} suffixed sequences → {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
