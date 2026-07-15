#!/usr/bin/env python3
"""
pipeline_summary.py — summarize sequence counts at each step of wvdb_build.

Reads key output FASTAs and TSVs, counts sequences, and writes:
  pipeline_summary.tsv      — machine-readable counts table (pandas DataFrame)
  pipeline_summary.md       — human-readable markdown report

Usage:
  pipeline_summary.py \\
      --all-fasta             filtered_all.fasta \\
      --clusters-tsv          vclust_clusters.tsv \\
      --trimming-candidates   trimming_candidates.fasta \\
      --trim-log              trim_selection.tsv \\
      --cluster-trim-complete complete_reps.fasta \\
      --blast-trim-input      blast_trim_input.fasta \\
      --blast-trim-complete   cluster_reps_complete.fasta \\
      --unvalidated-fasta     unvalidated_genomes.fasta \\
      --unvalidated-report    unvalidated_report.tsv \\
      --recluster-input       recluster_input.fasta \\
      --centroids             vclust_centroids.fasta \\
      --out-tsv               pipeline_summary.tsv \\
      --out-md                pipeline_summary.md \\
      --run-date              "2026-07-01" \\
      --ani                   0.95 \\
      --completeness          90
"""

import argparse
import sys
from datetime import date
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def count_fasta(path):
    """Return (n_seqs, total_bp, min_len, mean_len, max_len) or None if missing/empty."""
    if not path:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    lengths = [len(r.seq) for r in SeqIO.parse(path, 'fasta')]
    if not lengths:
        return None
    return (
        len(lengths),
        sum(lengths),
        min(lengths),
        round(sum(lengths) / len(lengths)),
        max(lengths),
    )


def count_clusters(tsv_path):
    """Return (n_sequences, n_clusters, n_singletons, n_multi) from vclust TSV."""
    p = Path(tsv_path)
    if not p.exists() or p.stat().st_size == 0:
        return (0, 0, 0, 0)
    df = pd.read_csv(tsv_path, sep='\t', comment='#')
    if df.shape[1] < 2:
        return (0, 0, 0, 0)
    cluster_col = df.columns[1]
    sizes = df[cluster_col].value_counts()
    n_clusters   = len(sizes)
    n_singletons = int((sizes == 1).sum())
    n_multi      = int((sizes > 1).sum())
    return (len(df), n_clusters, n_singletons, n_multi)


def count_unvalidated(report_path):
    """Return {reason: count} from unvalidated_report.tsv."""
    p = Path(report_path) if report_path else None
    if not p or not p.exists() or p.stat().st_size == 0:
        return {}
    df = pd.read_csv(report_path, sep='\t')
    if 'reason' not in df.columns:
        return {}
    return df['reason'].value_counts().to_dict()


def load_trim_log(path):
    """Return {source: count} from trim_selection.tsv."""
    defaults = {'trim12': 0, 'trim13': 0, 'trim23': 0}
    if not path:
        return defaults
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return defaults
    df = pd.read_csv(p, sep='\t')
    if 'source' not in df.columns:
        return defaults
    counts = df['source'].value_counts().to_dict()
    defaults.update(counts)
    return defaults


def fasta_row(step, description, stats):
    """Build a dict row from count_fasta() result tuple."""
    if stats is None:
        return dict(step=step, description=description,
                    n_seqs=0, total_bp=0, min_len=0, mean_len=0, max_len=0)
    return dict(step=step, description=description,
                n_seqs=stats[0], total_bp=stats[1],
                min_len=stats[2], mean_len=stats[3], max_len=stats[4])


def count_row(step, description, n_seqs):
    """Build a dict row for count-only entries (no length stats)."""
    return dict(step=step, description=description,
                n_seqs=n_seqs, total_bp=None, min_len=None,
                mean_len=None, max_len=None)


def fmt(n):
    """Format integer with thousands separator; None → '—'."""
    if n is None:
        return '—'
    return f"{int(n):,}"


def pct(part, total):
    if not total:
        return 'n/a'
    return f"{100 * part / total:.1f}%"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Summarize wvdb_build sequence counts at each pipeline step."
    )
    p.add_argument('--all-fasta',             required=True)
    p.add_argument('--clusters-tsv',          required=True)
    p.add_argument('--trimming-candidates',   default=None)
    p.add_argument('--trim-log',              default=None,
                   help='trim_selection.tsv from PICK_BEST_TRIM')
    p.add_argument('--cluster-trim-complete', default=None)
    p.add_argument('--blast-trim-input',      default=None)
    p.add_argument('--blast-trimmed',         default=None)
    p.add_argument('--blast-trim-complete',   default=None)
    p.add_argument('--unvalidated-fasta',     default=None)
    p.add_argument('--unvalidated-report',    default=None)
    p.add_argument('--recluster-input',       default=None)
    p.add_argument('--centroids',             required=True)
    p.add_argument('--out-tsv',               required=True)
    p.add_argument('--out-md',                required=True)
    p.add_argument('--run-date',              default=str(date.today()))
    p.add_argument('--ani',                   default='0.95')
    p.add_argument('--completeness',          default='90')
    args = p.parse_args()

    # --- Collect counts ---
    all_stats      = count_fasta(args.all_fasta)
    cluster_counts = count_clusters(args.clusters_tsv)
    cand_stats     = count_fasta(args.trimming_candidates)
    trim_src       = load_trim_log(args.trim_log)
    ct_complete    = count_fasta(args.cluster_trim_complete)
    bt_input       = count_fasta(args.blast_trim_input)
    bt_complete    = count_fasta(args.blast_trim_complete)
    unval_stats    = count_fasta(args.unvalidated_fasta)
    unval_reasons  = count_unvalidated(args.unvalidated_report)
    reclust_stats  = count_fasta(args.recluster_input)
    centroid_stats = count_fasta(args.centroids)

    n_input        = all_stats[0]      if all_stats      else 0
    n_clusters     = cluster_counts[1]
    n_singletons   = cluster_counts[2]
    n_multi        = cluster_counts[3]
    n_ct_complete  = ct_complete[0]    if ct_complete    else 0
    n_bt_input     = bt_input[0]       if bt_input       else 0
    n_bt_complete  = bt_complete[0]    if bt_complete    else 0
    n_unval        = unval_stats[0]    if unval_stats    else 0
    n_reclust      = reclust_stats[0]  if reclust_stats  else 0
    n_centroids    = centroid_stats[0] if centroid_stats else 0
    unval_no_hit   = unval_reasons.get('no_blast_hit', 0)
    unval_incomp   = unval_reasons.get('incomplete_after_trimming', 0)

    # --- Build DataFrame ---
    rows = [
        fasta_row('1_collect',              'All input genomes (virus + provirus)',
                  all_stats),
        count_row('2_clusters',             'vclust clusters',
                  n_clusters),
        count_row('2_singletons',           '  Singletons (1-member clusters)',
                  n_singletons),
        count_row('2_nonsingletons',        '  Non-singleton clusters',
                  n_multi),
        fasta_row('3_trimming_candidates',  'Rank1/2/3 trimming candidates',
                  cand_stats),
        count_row('3_trim_source_trim12',   '  Best trimming source: trim12',
                  trim_src['trim12']),
        count_row('3_trim_source_trim13',   '  Best trimming source: trim13',
                  trim_src['trim13']),
        count_row('3_trim_source_trim23',   '  Best trimming source: trim23',
                  trim_src['trim23']),
        fasta_row('3_cluster_trim_complete','Cluster-trimmed complete representatives',
                  ct_complete),
        fasta_row('4_blast_trim_input',     'Singletons + cluster-trim incomplete → blast_trim',
                  bt_input),
        fasta_row('4_blast_trim_complete',  'BLAST-trimmed complete representatives',
                  bt_complete),
        count_row('unvalidated_no_hit',     'Unvalidated: no BLAST hit',
                  unval_no_hit),
        count_row('unvalidated_incomplete', 'Unvalidated: incomplete after trimming',
                  unval_incomp),
        fasta_row('unvalidated_total',      'Unvalidated total',
                  unval_stats),
        fasta_row('5_recluster_input',      'All complete reps entering reclustering',
                  reclust_stats),
        fasta_row('5_centroids',            'Final non-redundant centroids',
                  centroid_stats),
    ]

    df = pd.DataFrame(rows, columns=[
        'step', 'description', 'n_seqs', 'total_bp', 'min_len', 'mean_len', 'max_len'
    ])
    df.to_csv(args.out_tsv, sep='\t', index=False)

    # --- Write markdown report ---
    md = f"""# wvdb_build Pipeline Summary

**Run date:** {args.run_date}
**ANI threshold:** {args.ani}  |  **Completeness threshold:** {args.completeness}%

---

## Sequence counts by step

| Step | Description | Sequences | Total bp | Mean len |
|------|-------------|----------:|----------:|---------:|
| 1. Collect | All input genomes (virus + provirus) | {fmt(n_input)} | {fmt(all_stats[1] if all_stats else None)} | {fmt(all_stats[3] if all_stats else None)} bp |
| 2. Cluster | vclust clusters formed | {fmt(n_clusters)} | — | — |
| | &nbsp;&nbsp;└ non-singleton clusters | {fmt(n_multi)} | — | — |
| | &nbsp;&nbsp;└ singletons | {fmt(n_singletons)} | — | — |
| 3. Cluster trim | Rank1/2/3 trimming candidates | {fmt(cand_stats[0] if cand_stats else None)} | — | — |
| | &nbsp;&nbsp;└ best source: trim12 | {fmt(trim_src['trim12'])} | — | — |
| | &nbsp;&nbsp;└ best source: trim13 | {fmt(trim_src['trim13'])} | — | — |
| | &nbsp;&nbsp;└ best source: trim23 | {fmt(trim_src['trim23'])} | — | — |
| | Complete trimmed reps | {fmt(n_ct_complete)} | {fmt(ct_complete[1] if ct_complete else None)} | {fmt(ct_complete[3] if ct_complete else None)} bp |
| 4. BLAST trim | Input (singletons + incomplete) | {fmt(n_bt_input)} | — | — |
| | Complete trimmed reps | {fmt(n_bt_complete)} | {fmt(bt_complete[1] if bt_complete else None)} | {fmt(bt_complete[3] if bt_complete else None)} bp |
| Unvalidated | No BLAST hit | {fmt(unval_no_hit)} | — | — |
| | Incomplete after trimming | {fmt(unval_incomp)} | — | — |
| | Total unvalidated | {fmt(n_unval)} | {fmt(unval_stats[1] if unval_stats else None)} | {fmt(unval_stats[3] if unval_stats else None)} bp |
| 5. Recluster | Input (cluster + blast trim complete) | {fmt(n_reclust)} | {fmt(reclust_stats[1] if reclust_stats else None)} | {fmt(reclust_stats[3] if reclust_stats else None)} bp |
| **5. Final** | **Non-redundant centroids** | **{fmt(n_centroids)}** | **{fmt(centroid_stats[1] if centroid_stats else None)}** | **{fmt(centroid_stats[3] if centroid_stats else None)} bp** |

---

## Retention summary

| Metric | Value |
|--------|------:|
| Input sequences | {fmt(n_input)} |
| Validated complete reps (cluster + blast trim) | {fmt(n_ct_complete + n_bt_complete)} ({pct(n_ct_complete + n_bt_complete, n_input)}) |
| Unvalidated sequences | {fmt(n_unval)} ({pct(n_unval, n_input)}) |
| Final centroids after reclustering | {fmt(n_centroids)} ({pct(n_centroids, n_input)}) |
| Redundancy removed by reclustering | {fmt((n_ct_complete + n_bt_complete) - n_centroids)} sequences |

---

*Generated by wvdb_build pipeline_summary.py*
"""

    with open(args.out_md, 'w') as f:
        f.write(md)

    print(f"[pipeline_summary] Written: {args.out_tsv}", file=sys.stderr)
    print(f"[pipeline_summary] Written: {args.out_md}", file=sys.stderr)


if __name__ == '__main__':
    main()
