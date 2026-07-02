#!/usr/bin/env python3
"""
pipeline_summary.py — summarize sequence counts at each step of wvdb_build.

Reads key output FASTAs and TSVs, counts sequences, and writes:
  pipeline_summary.tsv      — machine-readable counts table
  pipeline_summary.md       — human-readable markdown report

Usage:
  pipeline_summary.py \\
      --all-fasta             filtered_all.fasta \\
      --clusters-tsv          vclust_clusters.tsv \\
      --trimming-candidates   trimming_candidates.fasta \\
      --cluster-trim-complete complete_reps.fasta \\
      --blast-trim-input      blast_trim_input.fasta \\
      --blast-trimmed         blast_trimmed.fasta \\
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
    """Return (n_seqs, total_bp, min_len, mean_len, max_len) or None if missing."""
    if not path:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return (0, 0, 0, 0, 0)
    lengths = [len(r.seq) for r in SeqIO.parse(path, 'fasta')]
    if not lengths:
        return (0, 0, 0, 0, 0)
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
    # vclust format: object <TAB> cluster
    if df.shape[1] < 2:
        return (0, 0, 0, 0)
    cluster_col = df.columns[1]
    sizes = df[cluster_col].value_counts()
    n_clusters   = len(sizes)
    n_singletons = (sizes == 1).sum()
    n_multi      = (sizes > 1).sum()
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


def fmt(n):
    """Format integer with thousands separator."""
    if n is None:
        return 'n/a'
    return f"{n:,}"


def pct(part, total):
    """Format as percentage string."""
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
    ct_complete    = count_fasta(args.cluster_trim_complete)
    bt_input       = count_fasta(args.blast_trim_input)
    bt_trimmed     = count_fasta(args.blast_trimmed)
    bt_complete    = count_fasta(args.blast_trim_complete)
    unval_stats    = count_fasta(args.unvalidated_fasta)
    unval_reasons  = count_unvalidated(args.unvalidated_report)
    reclust_stats  = count_fasta(args.recluster_input)
    centroid_stats = count_fasta(args.centroids)

    n_input        = all_stats[0]    if all_stats    else 0
    n_clusters     = cluster_counts[1]
    n_singletons   = cluster_counts[2]
    n_multi        = cluster_counts[3]
    n_cand         = cand_stats[0]   if cand_stats   else 0
    n_ct_complete  = ct_complete[0]  if ct_complete  else 0
    n_bt_input     = bt_input[0]     if bt_input     else 0
    n_bt_complete  = bt_complete[0]  if bt_complete  else 0
    n_unval        = unval_stats[0]  if unval_stats  else 0
    n_reclust      = reclust_stats[0] if reclust_stats else 0
    n_centroids    = centroid_stats[0] if centroid_stats else 0

    # --- Write TSV ---
    rows = [
        ('step', 'description', 'n_seqs', 'total_bp', 'min_len', 'mean_len', 'max_len'),
        ('1_collect',            'All input genomes (virus + provirus)',
         *all_stats) if all_stats else ('1_collect', 'All input genomes', 0,0,0,0,0),
        ('2_cluster_candidates', 'Rank1/2/3 sequences sent to minimap2',
         *cand_stats) if cand_stats else ('2_cluster_candidates','',0,0,0,0,0),
        ('3_cluster_trim_complete', 'Cluster-trimmed complete representatives',
         *ct_complete) if ct_complete else ('3_cluster_trim_complete','',0,0,0,0,0),
        ('4_blast_trim_input',   'Singletons + cluster-trim incomplete → blast_trim',
         *bt_input) if bt_input else ('4_blast_trim_input','',0,0,0,0,0),
        ('4_blast_trim_complete','BLAST-trimmed complete representatives',
         *bt_complete) if bt_complete else ('4_blast_trim_complete','',0,0,0,0,0),
        ('unvalidated',          'No BLAST hit or incomplete after trimming',
         *unval_stats) if unval_stats else ('unvalidated','',0,0,0,0,0),
        ('5_recluster_input',    'All complete reps entering reclustering',
         *reclust_stats) if reclust_stats else ('5_recluster_input','',0,0,0,0,0),
        ('5_centroids',          'Final non-redundant centroids',
         *centroid_stats) if centroid_stats else ('5_centroids','',0,0,0,0,0),
    ]

    with open(args.out_tsv, 'w') as f:
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')

    # --- Write markdown report ---
    unval_no_hit    = unval_reasons.get('no_blast_hit', 0)
    unval_incomp    = unval_reasons.get('incomplete_after_trimming', 0)

    md = f"""# wvdb_build Pipeline Summary

**Run date:** {args.run_date}
**ANI threshold:** {args.ani}  |  **Completeness threshold:** {args.completeness}%

---

## Sequence counts by step

| Step | Description | Sequences | Total bp | Mean len |
|------|-------------|----------:|----------:|---------:|
| 1. Collect | All input genomes (virus + provirus) | {fmt(n_input)} | {fmt(all_stats[1] if all_stats else 0)} | {fmt(all_stats[3] if all_stats else 0)} bp |
| 2. Cluster | vclust clusters formed | {fmt(n_clusters)} clusters | — | — |
| 2. Cluster | &nbsp;&nbsp;└ multi-member clusters | {fmt(n_multi)} | — | — |
| 2. Cluster | &nbsp;&nbsp;└ singletons | {fmt(n_singletons)} | — | — |
| 3. Cluster trim | Rank1/2/3 candidates (minimap2 input) | {fmt(n_cand)} | — | — |
| 3. Cluster trim | Complete trimmed reps | {fmt(n_ct_complete)} | {fmt(ct_complete[1] if ct_complete else 0)} | {fmt(ct_complete[3] if ct_complete else 0)} bp |
| 4. BLAST trim | Input (singletons + incomplete) | {fmt(n_bt_input)} | — | — |
| 4. BLAST trim | Complete trimmed reps | {fmt(n_bt_complete)} | {fmt(bt_complete[1] if bt_complete else 0)} | {fmt(bt_complete[3] if bt_complete else 0)} bp |
| Unvalidated | No BLAST hit | {fmt(unval_no_hit)} | — | — |
| Unvalidated | Incomplete after trimming | {fmt(unval_incomp)} | — | — |
| 5. Recluster | Input (cluster_trim + blast_trim complete) | {fmt(n_reclust)} | {fmt(reclust_stats[1] if reclust_stats else 0)} | {fmt(reclust_stats[3] if reclust_stats else 0)} bp |
| **5. Final** | **Non-redundant centroids** | **{fmt(n_centroids)}** | **{fmt(centroid_stats[1] if centroid_stats else 0)}** | **{fmt(centroid_stats[3] if centroid_stats else 0)} bp** |

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
