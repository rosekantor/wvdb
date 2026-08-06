#!/usr/bin/env python3
"""
collect_unvalidated.py — collect all sequences that could not be validated
as high-quality trimmed complete genomes into a single unvalidated FASTA
with a summary report.

Sources of unvalidated sequences:
  --no-hit-ids       query IDs with no BLAST hit in any database
  --incomplete-ids   query IDs that were trimmed but CheckV incomplete
  --all-fasta        original untrimmed FASTA to retrieve sequences from

All sequences are written in their UNTRIMMED form — these are sequences for
which we cannot confirm chimera removal.

Usage:
  collect_unvalidated.py \\
      --no-hit-ids     no_hit_ids.txt \\
      --incomplete-ids blast_trim_incomplete_ids.txt \\
      --all-fasta      filtered_all.fasta \\
      --out-fasta      unvalidated_genomes.fasta \\
      --out-report     unvalidated_report.tsv
"""

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
import pandas as pd


def load_ids(path):
    """Load IDs from a file, one per line. Return empty set if missing."""
    if not path:
        return set()
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return set()
    with open(p) as f:
        return {line.strip() for line in f if line.strip()}


def main():
    parser = argparse.ArgumentParser(
        description="Collect unvalidated sequences into a single FASTA with report."
    )
    parser.add_argument('--no-hit-ids',     default=None,
                        help='IDs with no BLAST hit in any db')
    parser.add_argument('--incomplete-ids', default=None,
                        help='IDs trimmed but incomplete after CheckV')
    parser.add_argument('--all-fasta',      required=True,
                        help='Original untrimmed FASTA (filtered_all.fasta)')
    parser.add_argument('--out-fasta',      required=True,
                        help='Output FASTA of unvalidated sequences')
    parser.add_argument('--out-report',     required=True,
                        help='Output TSV report: seq_id, reason')
    args = parser.parse_args()

    no_hit_ids     = load_ids(args.no_hit_ids)
    incomplete_ids = load_ids(args.incomplete_ids)

    # Build reason map
    reason = {}
    for id_ in no_hit_ids:
        reason[id_] = 'no_qualifying_hit'
    for id_ in incomplete_ids:
        # incomplete takes priority if somehow in both sets
        reason[id_] = 'incomplete_after_trimming'

    all_unvalidated = set(reason.keys())

    if not all_unvalidated:
        print("[collect_unvalidated] no unvalidated sequences", file=sys.stderr)
        open(args.out_fasta,  'w').close()
        open(args.out_report, 'w').close()
        return

    # Retrieve untrimmed sequences from all_fasta
    written = 0
    not_found = []
    with open(args.out_fasta, 'w') as fout:
        for rec in SeqIO.parse(args.all_fasta, 'fasta'):
            if rec.id in all_unvalidated:
                SeqIO.write(rec, fout, 'fasta')
                written += 1

    not_found = all_unvalidated - {
        rec.id for rec in SeqIO.parse(args.out_fasta, 'fasta')
    }

    # Write report
    rows = [{'seq_id': id_, 'reason': r} for id_, r in sorted(reason.items())]
    pd.DataFrame(rows).to_csv(args.out_report, sep='\t', index=False)

    print(f"[collect_unvalidated] unvalidated sequences  : {len(all_unvalidated)}",
          file=sys.stderr)
    print(f"[collect_unvalidated]   no qualifying hit    : {len(no_hit_ids)}",
          file=sys.stderr)
    print(f"[collect_unvalidated]   incomplete after trim: {len(incomplete_ids)}",
          file=sys.stderr)
    print(f"[collect_unvalidated] written to FASTA       : {written}",
          file=sys.stderr)
    if not_found:
        print(f"[collect_unvalidated] WARNING: {len(not_found)} IDs not found "
              f"in all_fasta", file=sys.stderr)


if __name__ == '__main__':
    main()
