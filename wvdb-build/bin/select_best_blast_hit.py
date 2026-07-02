#!/usr/bin/env python3
"""
select_best_blast_hit.py — select the best BLAST reference per query across
initial and secondary databases, preferring the initial db when both have hits.

Reads ANI TSVs from both databases and outputs three files:
  initial_hits.tsv     — queries assigned to initial db (initial had a hit)
  secondary_hits.tsv   — queries assigned to secondary db (initial had no hit,
                          secondary did)
  no_hit_ids.txt       — queries with no hit in either db → unvalidated

Usage:
  select_best_blast_hit.py \\
      --initial-ani   initial.blastn.ani.tsv \\
      --secondary-ani secondary.blastn.ani.tsv  (optional) \\
      --initial-out   initial_hits.tsv \\
      --secondary-out secondary_hits.tsv \\
      --no-hit-ids    no_hit_ids.txt \\
      --all-query-ids blast_trim_input_ids.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_ani(path):
    """Load ANI TSV, return empty DataFrame if missing or empty."""
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(p, sep='\t')


def load_ids(path):
    """Load one-per-line IDs file into a set."""
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return set()
    with open(p) as f:
        return {line.strip() for line in f if line.strip()}


def main():
    parser = argparse.ArgumentParser(
        description="Route blast_trim queries to initial or secondary db hit."
    )
    parser.add_argument('--initial-ani',   required=True,
                        help='ANI TSV from blastani for initial db')
    parser.add_argument('--secondary-ani', default=None,
                        help='ANI TSV from blastani for secondary db (optional)')
    parser.add_argument('--initial-out',   required=True,
                        help='Output ANI TSV for initial db queries')
    parser.add_argument('--secondary-out', required=True,
                        help='Output ANI TSV for secondary db queries')
    parser.add_argument('--no-hit-ids',    required=True,
                        help='Output file: query IDs with no hit in either db')
    parser.add_argument('--all-query-ids', required=True,
                        help='File listing all query IDs (from blast_trim_input.fasta)')
    args = parser.parse_args()

    all_queries = load_ids(args.all_query_ids)
    initial_ani = load_ani(args.initial_ani)
    secondary_ani = load_ani(args.secondary_ani)

    # Get query IDs with a hit in each db
    initial_hits = set()
    if not initial_ani.empty:
        q_col = initial_ani.columns[0]
        initial_hits = set(initial_ani[q_col].astype(str))

    secondary_hits = set()
    if not secondary_ani.empty:
        q_col = secondary_ani.columns[0]
        secondary_hits = set(secondary_ani[q_col].astype(str))

    # Route: initial preferred; secondary only if no initial hit
    assigned_initial   = initial_hits
    assigned_secondary = secondary_hits - initial_hits
    no_hit             = all_queries - initial_hits - secondary_hits

    # Write initial hits TSV (full rows)
    if not initial_ani.empty:
        q_col = initial_ani.columns[0]
        initial_ani[initial_ani[q_col].isin(assigned_initial)] \
            .to_csv(args.initial_out, sep='\t', index=False)
    else:
        open(args.initial_out, 'w').close()

    # Write secondary hits TSV (only queries not in initial)
    if not secondary_ani.empty and assigned_secondary:
        q_col = secondary_ani.columns[0]
        secondary_ani[secondary_ani[q_col].isin(assigned_secondary)] \
            .to_csv(args.secondary_out, sep='\t', index=False)
    else:
        open(args.secondary_out, 'w').close()

    # Write no-hit IDs
    with open(args.no_hit_ids, 'w') as f:
        for id_ in sorted(no_hit):
            f.write(id_ + '\n')

    print(f"[select_best_blast_hit] total queries     : {len(all_queries)}",
          file=sys.stderr)
    print(f"[select_best_blast_hit] → initial db      : {len(assigned_initial)}",
          file=sys.stderr)
    print(f"[select_best_blast_hit] → secondary db    : {len(assigned_secondary)}",
          file=sys.stderr)
    print(f"[select_best_blast_hit] → no hit          : {len(no_hit)} → unvalidated",
          file=sys.stderr)


if __name__ == '__main__':
    main()
