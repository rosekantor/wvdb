#!/usr/bin/env python3
"""
select_best_blast_hit.py — select the best BLAST reference per query across
initial and secondary databases, preferring the initial db when both have
a qualifying hit.

A hit "qualifies" for trimming if it meets BOTH:
  target coverage (tcov) >= --min-tcov   (default: 85%)
  percent identity (pid) >= --min-pid    (default: 85%, genus-level)

Target coverage (not query coverage) is used because the untrimmed query
may be larger than the true genome due to chimeric/host-contaminant
regions — full query coverage of the hit is not expected. The identity
threshold ensures trimming only happens against a genuinely homologous
reference, not a spurious low-identity hit.

Queries with a hit that exists but does not clear these thresholds are
treated the same as queries with no hit at all — both are combined into
a single "no_qualifying_hit" category and routed to unvalidated. For
this dataset (novel viruses), the vast majority of "hits" are very short,
low-coverage matches that are not meaningfully different from no hit.

Reads ANI TSVs from both databases and outputs three files:
  initial_hits.tsv     — queries with a qualifying hit in the initial db
  secondary_hits.tsv   — queries with a qualifying hit in the secondary db
                          only (no qualifying initial hit)
  no_hit_ids.txt       — queries with no qualifying hit in either db
                          → unvalidated (reason: no_qualifying_hit)

Usage:
  select_best_blast_hit.py \\
      --initial-ani   initial.blastn.ani.tsv \\
      --secondary-ani secondary.blastn.ani.tsv  (optional) \\
      --initial-out   initial_hits.tsv \\
      --secondary-out secondary_hits.tsv \\
      --no-hit-ids    no_hit_ids.txt \\
      --all-query-ids blast_trim_input_ids.txt \\
      --min-tcov      85 \\
      --min-pid       85
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


def filter_qualifying(ani_df, min_tcov, min_pid):
    """
    Restrict an ANI DataFrame to rows meeting the qualifying-hit thresholds:
    tcov >= min_tcov AND pid >= min_pid. Returns an empty DataFrame with the
    same columns if input is empty or no rows qualify.
    """
    if ani_df.empty:
        return ani_df
    missing = {'tcov', 'pid'} - set(ani_df.columns)
    if missing:
        print(f"ERROR: ANI TSV missing required column(s): {missing}",
              file=sys.stderr)
        sys.exit(1)
    return ani_df[(ani_df.tcov >= min_tcov) & (ani_df.pid >= min_pid)].copy()


def main():
    parser = argparse.ArgumentParser(
        description="Route blast_trim queries to initial or secondary db hit, "
                    "requiring a qualifying hit (tcov + pid thresholds)."
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
                        help='Output file: query IDs with no qualifying hit in either db')
    parser.add_argument('--all-query-ids', required=True,
                        help='File listing all query IDs (from blast_trim_input.fasta)')
    parser.add_argument('--min-tcov', type=float, default=85.0,
                        help='Minimum target coverage %% for a hit to qualify (default: 85)')
    parser.add_argument('--min-pid',  type=float, default=85.0,
                        help='Minimum percent identity for a hit to qualify (default: 85, genus-level)')
    args = parser.parse_args()

    all_queries   = load_ids(args.all_query_ids)
    initial_ani   = load_ani(args.initial_ani)
    secondary_ani = load_ani(args.secondary_ani)

    n_initial_raw   = len(initial_ani)
    n_secondary_raw = len(secondary_ani)

    initial_qual   = filter_qualifying(initial_ani,   args.min_tcov, args.min_pid)
    secondary_qual = filter_qualifying(secondary_ani, args.min_tcov, args.min_pid)

    # Get query IDs with a QUALIFYING hit in each db
    initial_hits = set()
    if not initial_qual.empty:
        q_col = initial_qual.columns[0]
        initial_hits = set(initial_qual[q_col].astype(str))

    secondary_hits = set()
    if not secondary_qual.empty:
        q_col = secondary_qual.columns[0]
        secondary_hits = set(secondary_qual[q_col].astype(str))

    # Route: initial preferred; secondary only if no qualifying initial hit
    assigned_initial   = initial_hits
    assigned_secondary = secondary_hits - initial_hits
    no_hit              = all_queries - initial_hits - secondary_hits

    # Write initial hits TSV (full rows, qualifying only)
    if not initial_qual.empty:
        q_col = initial_qual.columns[0]
        initial_qual[initial_qual[q_col].isin(assigned_initial)] \
            .to_csv(args.initial_out, sep='\t', index=False)
    else:
        open(args.initial_out, 'w').close()

    # Write secondary hits TSV (only queries not qualifying in initial)
    if not secondary_qual.empty and assigned_secondary:
        q_col = secondary_qual.columns[0]
        secondary_qual[secondary_qual[q_col].isin(assigned_secondary)] \
            .to_csv(args.secondary_out, sep='\t', index=False)
    else:
        open(args.secondary_out, 'w').close()

    # Write no-hit IDs (true no-hit + hit-but-not-qualifying, combined)
    with open(args.no_hit_ids, 'w') as f:
        for id_ in sorted(no_hit):
            f.write(id_ + '\n')

    print(f"[select_best_blast_hit] total queries          : {len(all_queries)}",
          file=sys.stderr)
    print(f"[select_best_blast_hit] initial db  : {n_initial_raw} raw hits → "
          f"{len(initial_qual)} qualifying (tcov>={args.min_tcov}, pid>={args.min_pid})",
          file=sys.stderr)
    print(f"[select_best_blast_hit] secondary db : {n_secondary_raw} raw hits → "
          f"{len(secondary_qual)} qualifying",
          file=sys.stderr)
    print(f"[select_best_blast_hit] → initial db      : {len(assigned_initial)}",
          file=sys.stderr)
    print(f"[select_best_blast_hit] → secondary db    : {len(assigned_secondary)}",
          file=sys.stderr)
    print(f"[select_best_blast_hit] → no_qualifying_hit: {len(no_hit)} → unvalidated",
          file=sys.stderr)


if __name__ == '__main__':
    main()
