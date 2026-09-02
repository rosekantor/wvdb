#!/usr/bin/env python3
"""
annotation_summary.py — summarize annotation coverage across the wvdb_annotate
pipeline outputs.

Reads CheckV, geNomad, the merged annotation table, and (optionally) host
predictions, and writes:
  annotation_summary.tsv   — machine-readable counts table (pandas DataFrame)
  annotation_summary.md    — human-readable markdown report

Usage:
  annotation_summary.py \\
      --fasta            votus.fasta \\
      --checkv           checkv_out/quality_summary.tsv \\
      --genomad          genomad_out/<prefix>_virus_summary.tsv \\
      --merged           merged_annotations.tsv \\
      [--host-predictions host_predictions.tsv] \\
      --out-tsv          annotation_summary.tsv \\
      --out-md           annotation_summary.md \\
      --run-date         "2026-07-30"
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

def count_fasta_seqs(path):
    """Return sequence count in a FASTA file, or None if missing/empty."""
    if not path:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    return sum(1 for _ in SeqIO.parse(path, 'fasta'))


def load_tsv(path, required_cols=None):
    """Load a TSV, return None if missing/empty. Warn if required cols absent."""
    if not path:
        return None
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    df = pd.read_csv(path, sep='\t')
    if required_cols:
        missing = set(required_cols) - set(df.columns)
        if missing:
            print(f"Warning: {path} missing expected columns: {missing}",
                  file=sys.stderr)
    return df


def count_row(step, description, n_seqs, pct_of=None, total=None):
    """Build a dict row with an optional percentage-of-total annotation."""
    row = dict(step=step, description=description, n_seqs=n_seqs)
    if pct_of is not None and total:
        row['pct_total'] = round(100 * pct_of / total, 1)
    else:
        row['pct_total'] = None
    return row


def fmt(n):
    if n is None:
        return '—'
    return f"{int(n):,}"


def pct(part, total):
    if not total:
        return 'n/a'
    return f"{100 * part / total:.1f}%"


def find_blastn_hit_cols(df):
    """Return sorted list of database names from hit_level_<db> columns."""
    cols = [c for c in df.columns if c.startswith('hit_level_')]
    return sorted(c.replace('hit_level_', '') for c in cols)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Summarize wvdb_annotate annotation coverage."
    )
    p.add_argument('--fasta',            required=True)
    p.add_argument('--checkv',           required=True)
    p.add_argument('--genomad',          required=True)
    p.add_argument('--merged',           required=True)
    p.add_argument('--host-predictions', default=None)
    p.add_argument('--out-tsv',          required=True)
    p.add_argument('--out-md',           required=True)
    p.add_argument('--run-date',         default=str(date.today()))
    args = p.parse_args()

    n_input = count_fasta_seqs(args.fasta) or 0

    checkv_df  = load_tsv(args.checkv,  required_cols=['contig_id', 'checkv_quality'])
    genomad_df = load_tsv(args.genomad, required_cols=['seq_name', 'taxonomy'])
    merged_df  = load_tsv(args.merged)
    host_df    = load_tsv(args.host_predictions)

    if merged_df is None or merged_df.empty:
        print("ERROR: merged_annotations.tsv is missing or empty", file=sys.stderr)
        sys.exit(1)

    n_merged = len(merged_df)

    # --- CheckV quality tiers ---
    checkv_tiers = {}
    if 'checkv_quality' in merged_df.columns:
        checkv_tiers = merged_df['checkv_quality'].value_counts().to_dict()

    # --- geNomad / RdRp consensus classification coverage ---
    n_order_classified  = 0
    n_family_classified = 0
    if 'Order_consensus' in merged_df.columns:
        n_order_classified = (merged_df['Order_consensus'] != 'Unclassified').sum()
    if 'Family_consensus' in merged_df.columns:
        n_family_classified = (merged_df['Family_consensus'] != 'Unclassified').sum()

    # --- RdRPCATCH coverage (Order_RdRp non-null) ---
    n_rdrp_hit = 0
    if 'Order_RdRp' in merged_df.columns:
        n_rdrp_hit = merged_df['Order_RdRp'].notna().sum()

    # --- BLASTn hit coverage per database ---
    db_names = find_blastn_hit_cols(merged_df)
    blastn_hits = {}
    for db in db_names:
        col = f'hit_level_{db}'
        blastn_hits[db] = (merged_df[col].isin(['genus', 'species'])).sum()

    # --- RNAVirHost coverage ---
    n_rnavirhost = 0
    if 'rnavirhost_domain' in merged_df.columns:
        n_rnavirhost = (merged_df['rnavirhost_domain'] != 'Unknown').sum()

    # --- ICTV host lookup coverage ---
    n_ictv_host = 0
    if 'host_ICTV_simple' in merged_df.columns:
        n_ictv_host = (merged_df['host_ICTV_simple'] != 'Unknown').sum()

    # --- Final host determination (if guess_host was run) ---
    n_final_host = 0
    host_breakdown = {}
    if host_df is not None and 'final_host_determination' in host_df.columns:
        n_final_host = (host_df['final_host_determination'] != 'unknown').sum()
        host_breakdown = host_df['final_host_determination'].value_counts().to_dict()

    # --- Build DataFrame ---
    rows = [
        count_row('1_input',            'Input vOTUs',                n_input),
        count_row('2_checkv_genomad',   'Sequences in merged annotation table',
                  n_merged, n_merged, n_input),
    ]
    for tier, n in sorted(checkv_tiers.items()):
        rows.append(count_row(f'2_checkv_{tier}', f'  CheckV quality: {tier}',
                              n, n, n_merged))

    rows.append(count_row('3_order_classified',  'Classified at Order level (RdRp/geNomad)',
                          n_order_classified, n_order_classified, n_merged))
    rows.append(count_row('3_family_classified', 'Classified at Family level (RdRp/geNomad)',
                          n_family_classified, n_family_classified, n_merged))
    rows.append(count_row('3_rdrp_hit',          'RdRPCATCH RdRp domain detected',
                          n_rdrp_hit, n_rdrp_hit, n_merged))

    for db in db_names:
        rows.append(count_row(f'4_blastn_{db}', f'  BLASTn hit (genus/species): {db}',
                              blastn_hits[db], blastn_hits[db], n_merged))

    rows.append(count_row('5_rnavirhost',    'RNAVirHost domain prediction (not Unknown)',
                          n_rnavirhost, n_rnavirhost, n_merged))
    rows.append(count_row('6_ictv_host',     'ICTV host lookup (Family_consensus matched)',
                          n_ictv_host, n_ictv_host, n_merged))

    if host_df is not None:
        rows.append(count_row('7_final_host', 'Final host determined (guess_host)',
                              n_final_host, n_final_host, n_merged))
        for host, n in sorted(host_breakdown.items(), key=lambda x: -x[1]):
            rows.append(count_row(f'7_host_{host}', f'  host_final: {host}',
                                  n, n, n_merged))

    df = pd.DataFrame(rows, columns=['step', 'description', 'n_seqs', 'pct_total'])
    df.to_csv(args.out_tsv, sep='\t', index=False)

    # --- Write markdown report ---
    checkv_lines = '\n'.join(
        f"| | &nbsp;&nbsp;└ {tier} | {fmt(n)} | {pct(n, n_merged)} |"
        for tier, n in sorted(checkv_tiers.items())
    )
    blastn_lines = '\n'.join(
        f"| 4. BLASTn | {db} (genus/species hit) | {fmt(blastn_hits[db])} | {pct(blastn_hits[db], n_merged)} |"
        for db in db_names
    ) if db_names else "| 4. BLASTn | (not run) | — | — |"

    host_section = ""
    if host_df is not None:
        host_lines = '\n'.join(
            f"| | &nbsp;&nbsp;└ {host} | {fmt(n)} | {pct(n, n_merged)} |"
            for host, n in sorted(host_breakdown.items(), key=lambda x: -x[1])
        )
        host_section = f"""
| 7. Final host | Determined (guess_host) | {fmt(n_final_host)} | {pct(n_final_host, n_merged)} |
{host_lines}
"""
    else:
        host_section = "\n| 7. Final host | (guess_host not run) | — | — |\n"

    md = f"""# wvdb_annotate Annotation Summary

**Run date:** {args.run_date}

---

## Annotation coverage by step

| Step | Description | Sequences | % of total |
|------|-------------|----------:|-----------:|
| 1. Input | Input vOTUs | {fmt(n_input)} | — |
| 2. CheckV/geNomad | In merged annotation table | {fmt(n_merged)} | {pct(n_merged, n_input)} |
{checkv_lines}
| 3. Classification | Order classified (RdRp/geNomad) | {fmt(n_order_classified)} | {pct(n_order_classified, n_merged)} |
| 3. Classification | Family classified (RdRp/geNomad) | {fmt(n_family_classified)} | {pct(n_family_classified, n_merged)} |
| 3. Classification | RdRPCATCH RdRp domain detected | {fmt(n_rdrp_hit)} | {pct(n_rdrp_hit, n_merged)} |
{blastn_lines}
| 5. RNAVirHost | Domain predicted (not Unknown) | {fmt(n_rnavirhost)} | {pct(n_rnavirhost, n_merged)} |
| 6. ICTV | Host lookup matched | {fmt(n_ictv_host)} | {pct(n_ictv_host, n_merged)} |
{host_section}
---

*Generated by wvdb_annotate annotation_summary.py*
"""

    with open(args.out_md, 'w') as f:
        f.write(md)

    print(f"[annotation_summary] Written: {args.out_tsv}", file=sys.stderr)
    print(f"[annotation_summary] Written: {args.out_md}", file=sys.stderr)


if __name__ == '__main__':
    main()
