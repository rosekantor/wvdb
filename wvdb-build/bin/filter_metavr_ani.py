#!/usr/bin/env python3
"""
filter_metavr_ani.py — remove from the metavr ANI TSV any queries that
were already trimmed by the refseq_ev database.

This ensures TRIM_GENOMES_BLAST_METAVR only runs nucmer on genuinely
metavr-only hits, avoiding redundant work on sequences already handled.

Usage:
  filter_metavr_ani.py \\
      --metavr-ani   metavr.blastn.ani.tsv \\
      --refseq-bed   refseq_ev.trimming.bed \\
      --out          metavr_ani_filtered.tsv
"""

import argparse
import sys
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--metavr-ani',  required=True, help='metavr ANI TSV from blastani')
    p.add_argument('--refseq-bed',  required=True, help='trimming.bed from refseq_ev trim step')
    p.add_argument('--out',         required=True, help='filtered output ANI TSV')
    args = p.parse_args()

    # Load query IDs already trimmed by refseq_ev from their BED file
    already_trimmed = set()
    with open(args.refseq_bed) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                already_trimmed.add(parts[0])

    # Load metavr ANI TSV and filter out already-trimmed queries
    ani = pd.read_csv(args.metavr_ani, sep='\t')
    if ani.empty:
        ani.to_csv(args.out, sep='\t', index=False)
        print("filter_metavr_ani: input ANI TSV is empty", file=sys.stderr)
        return

    query_col = ani.columns[0]   # first column is query name
    filtered = ani[~ani[query_col].isin(already_trimmed)]
    filtered.to_csv(args.out, sep='\t', index=False)

    n_removed = len(ani) - len(filtered)
    print(
        f"filter_metavr_ani: {len(already_trimmed)} refseq_ev trimmed IDs, "
        f"{n_removed} removed from metavr ANI, "
        f"{len(filtered)} metavr-only pairs remain",
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
