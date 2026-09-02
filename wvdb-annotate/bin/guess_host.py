#!/usr/bin/env python3
"""
guess_host.py — ensemble host determination combining ICTV taxonomy,
RNAVirHost prediction, and (optionally) categorized core_nt BLAST hit
metadata.

This script has no LLM dependency and is safe to rerun freely — the slow,
API-cost-incurring categorization of free-text host/isolation-source
terms happens separately in categorize_host_terms.py, whose cached output
(host-dict / isolation-dict TSVs) is consumed here.

Decision tree:
  * compare ICTV vs RNAVirHost, determine agreement
  * where BLAST-derived host exists:
      * ICTV-RNAVirHost agree, BLAST agrees too    → use BLAST host
      * ICTV-RNAVirHost agree, BLAST disagrees     → use RNAVirHost, flag for review
      * ICTV-RNAVirHost disagree, BLAST agrees ICTV → use BLAST host
      * ICTV-RNAVirHost disagree, no clear BLAST agreement → use ICTV, flag for review
      * ICTV-RNAVirHost both unknown               → use BLAST, flag if isolation
                                                       source is fecal (dietary risk)
  * where no BLAST-derived host exists:
      * ICTV-RNAVirHost agree    → use RNAVirHost
      * ICTV-RNAVirHost disagree → use whichever is known (ICTV preferred)

Manual corrections applied after the decision tree:
  * Unclassified Order + Class in a known-bacteriophage class
    (Leviviricetes, Caudoviricetes, Vidaverviricetes, Faserviricetes)
    → bacteria
  * Norzivirales / Timlovirales → always bacteria (well-established)
  * Martellivirales → flagged for review (usually but not always plants;
    written to manual_review.tsv, not auto-corrected)
  * Rare host categories (< --other-threshold vOTUs) consolidated to "other"

Usage:
  guess_host.py \\
      --merged              merged_annotations.tsv \\
      --host-dict           host_dict.tsv \\
      --isolation-dict      isolation_dict.tsv \\
      --out                 host_predictions.tsv \\
      --out-manual-review   manual_review.tsv \\
      --other-threshold     100
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# Bacteriophage classes where an Unclassified Order should still be called
# bacteria based on Class alone
BACTERIOPHAGE_CLASSES = ['Leviviricetes', 'Caudoviricetes', 'Vidaverviricetes', 'Faserviricetes']

# Orders that are essentially always one host type regardless of what the
# decision tree concluded — well-established exceptions, not guesses
ALWAYS_BACTERIA_ORDERS = ['Norzivirales', 'Timlovirales']

# Orders that are USUALLY one host type but not reliably enough to
# auto-correct — flagged for manual review instead
REVIEW_ORDERS = {
    'Martellivirales': 'plants',  # usually but not always plants
}


def load_dict_tsv(path):
    """Load a two-column TSV as {key: value}. Empty dict if missing."""
    if not path:
        return {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return {}
    df = pd.read_csv(p, sep='\t')
    if df.shape[1] < 2:
        return {}
    key_col, val_col = df.columns[0], df.columns[1]
    return dict(zip(df[key_col].astype(str), df[val_col]))


def ictv_vs_rnavirhost_agreement(row):
    """Compare RNAVirHost to ICTV. ICTV may contain multiple comma-separated values."""
    rna = str(row["rnavirhost_domain"]).strip()
    ictv = str(row["host_ICTV"]).strip()

    if rna == "unknown" or ictv == "unknown":
        return "unknown"

    ictv_vals = [x.strip() for x in ictv.split(",")]
    return "agree" if rna in ictv_vals else "disagree"


def out(host, tool, flag=False, reason=None):
    return pd.Series([host, tool, flag, reason])


def choose_final(row):
    rna = str(row["rnavirhost_domain"]).strip()
    ictv = str(row["host_ICTV"]).strip()
    ictv_vals = [] if ictv == "unknown" else [x.strip() for x in ictv.split(",")]

    blast_host = row["hosts_ntBlastHit_simple"]
    blast_iso = row["isolation_source_ntBlastHit_simple"] if pd.notna(row["isolation_source_ntBlastHit_simple"]) else None
    has_blast = pd.notna(blast_host)

    rna_known = rna != "unknown"
    ictv_known = ictv != "unknown"
    both_known = rna_known and ictv_known
    both_unknown = (not rna_known) and (not ictv_known)

    ictv_rna_agree = both_known and (rna in ictv_vals)
    blast_agrees_ictv = has_blast and ictv_known and (blast_host in ictv_vals)
    blast_agrees_rna = has_blast and rna_known and (blast_host == rna)

    if ictv_rna_agree:
        if blast_agrees_ictv:
            return out(blast_host, "ICTV_RNAVirHost_ntBlastHit_host")
        if has_blast and not blast_agrees_ictv:
            return out(rna, "ICTV_RNAVirHost", True, "ICTV and RNAVirHost agree, BLAST disagrees")
        return out(rna, "ICTV_RNAVirHost")

    if both_known and not ictv_rna_agree:
        if has_blast and blast_agrees_ictv:
            return out(blast_host, "ICTV_ntBlastHit_host")
        if not has_blast:
            return out(ictv, "ICTV")
        if blast_agrees_rna:
            return out(blast_host, "RNAVirHost_ntBlastHit_host")
        return out(ictv, "ICTV", True, "ICTV and RNAVirHost disagree, BLAST disagrees with both")

    if not both_known:
        if both_unknown and has_blast:
            if str(blast_iso) == "fecal associated":
                return out(blast_host, "ntBlastHit_host", True, "BLAST used alone, isolation source is fecal associated")
            return out(blast_host, "ntBlastHit_host")
        if both_unknown and not has_blast:
            return out("unknown", "none", False, None)
        if has_blast:
            if ictv_known and blast_agrees_ictv:
                return out(blast_host, "ICTV_ntBlastHit_host")
            if rna_known and blast_agrees_rna:
                return out(blast_host, "RNAVirHost_ntBlastHit_host")
            if ictv_known:
                return out(ictv, "ICTV")
            if not ictv_known:
                if str(blast_iso) == "fecal associated":
                    return out(blast_host, "ntBlastHit_host", True, "BLAST used alone, isolation source is fecal associated")
                return out(blast_host, "ntBlastHit_host")
        if not has_blast:
            if ictv_known:
                return out(ictv, "ICTV")
            if rna_known:
                return out(rna, "RNAVirHost")

    return out("review", "review", True, "Unhandled case")


def main():
    p = argparse.ArgumentParser(
        description="Ensemble host determination from ICTV, RNAVirHost, "
                    "and categorized core_nt BLAST hit metadata."
    )
    p.add_argument('--merged',            required=True,
                   help='merged_annotations.tsv from merge_annotations.py')
    p.add_argument('--host-dict',         default=None,
                   help='Categorized host dict TSV from categorize_host_terms.py')
    p.add_argument('--isolation-dict',    default=None,
                   help='Categorized isolation-source dict TSV from categorize_host_terms.py')
    p.add_argument('--out',               required=True)
    p.add_argument('--out-manual-review', required=True,
                   help='vOTUs needing manual review (e.g. Martellivirales), '
                        'kept separate from auto-corrected results')
    p.add_argument('--other-threshold',   type=int, default=100,
                   help='Host categories with fewer than this many vOTUs are '
                        'consolidated into "other" (default: 100)')
    args = p.parse_args()

    df = pd.read_csv(args.merged, sep='\t')

    if 'contig' not in df.columns:
        print("ERROR: merged_annotations.tsv missing 'contig' column", file=sys.stderr)
        sys.exit(1)

    df.loc[df['host_ICTV'].isna(), 'host_ICTV'] = 'unknown'
    df['rnavirhost_domain'] = df['rnavirhost_domain'].str.lower()
    df['host_ICTV'] = df['host_ICTV'].str.lower()

    # Normalize RNAVirHost's vocabulary to match ICTV's
    df.loc[df['rnavirhost_domain'] == 'viridiplantae',  'rnavirhost_domain'] = 'plants'
    df.loc[df['rnavirhost_domain'] == 'chordata',       'rnavirhost_domain'] = 'vertebrates'
    df.loc[df['rnavirhost_domain'] == 'invertebrate',   'rnavirhost_domain'] = 'invertebrates'

    # Merge in categorized BLAST host/isolation-source terms, if available
    host_dict      = load_dict_tsv(args.host_dict)
    isolation_dict = load_dict_tsv(args.isolation_dict)

    if 'hosts_ntBlastHit' in df.columns:
        df['hosts_ntBlastHit_simple'] = df['hosts_ntBlastHit'].map(host_dict)
    else:
        df['hosts_ntBlastHit_simple'] = np.nan

    if 'isolation_source_ntBlastHit' in df.columns:
        df['isolation_source_ntBlastHit_simple'] = df['isolation_source_ntBlastHit'].map(isolation_dict)
    else:
        df['isolation_source_ntBlastHit_simple'] = np.nan

    # Apply decision tree
    df["ictv_rnavirhost_agreement"] = df.apply(ictv_vs_rnavirhost_agreement, axis=1)
    df["blast_exists"] = df["hosts_ntBlastHit_simple"].notna()

    df[[
        "final_host_determination", "final_tool_used",
        "review_flag", "review_reason"
    ]] = df.apply(choose_final, axis=1)

    # Manual correction: Unclassified Order in a known bacteriophage Class
    unclassified_bacteriophage = (
        (df['Order_consensus'] == 'Unclassified') &
        (df.get('Class_gNd', pd.Series(dtype=str)).isin(BACTERIOPHAGE_CLASSES) |
         df.get('Class_RdRp', pd.Series(dtype=str)).isin(BACTERIOPHAGE_CLASSES))
    )
    df.loc[unclassified_bacteriophage, 'final_host_determination'] = 'bacteria'
    df.loc[unclassified_bacteriophage, 'final_tool_used'] = 'manual_by_class'

    # Manual correction: always-bacteria orders
    for order in ALWAYS_BACTERIA_ORDERS:
        mask = (df['Order_consensus'] == order) & (df['final_host_determination'] != 'bacteria')
        df.loc[mask, 'final_host_determination'] = 'bacteria'
        df.loc[mask, 'final_tool_used'] = 'manual_by_order'

    # Flag (don't auto-correct) review orders like Martellivirales
    manual_review_rows = []
    for order, expected_host in REVIEW_ORDERS.items():
        mask = (df['Order_consensus'] == order) & (df['final_host_determination'] != expected_host)
        if mask.any():
            review_subset = df.loc[mask, ['contig', 'Order_consensus', 'final_host_determination']].copy()
            review_subset['expected_host'] = expected_host
            review_subset['review_reason_manual'] = (
                f"Order {order} is usually {expected_host}-associated but "
                f"final determination differs — confirm manually"
            )
            manual_review_rows.append(review_subset)

    if manual_review_rows:
        manual_review_df = pd.concat(manual_review_rows, ignore_index=True)
    else:
        manual_review_df = pd.DataFrame(columns=[
            'contig', 'Order_consensus', 'final_host_determination',
            'expected_host', 'review_reason_manual'
        ])
    manual_review_df.to_csv(args.out_manual_review, sep='\t', index=False)

    # Consolidate rare host categories into "other"
    summary_host = df.groupby('final_host_determination')[['contig']].count().reset_index()
    rare = summary_host[summary_host['contig'] < args.other_threshold]['final_host_determination'].tolist()

    df['host_final'] = df['final_host_determination']
    df.loc[df['final_host_determination'].isin(rare), 'host_final'] = 'other'

    df.to_csv(args.out, sep='\t', index=False)

    n_flagged = df['review_flag'].sum()
    print(f"[guess_host] {len(df)} vOTUs processed", file=sys.stderr)
    print(f"[guess_host] {n_flagged} flagged for review (decision tree)", file=sys.stderr)
    print(f"[guess_host] {len(manual_review_df)} flagged for manual review "
          f"(ambiguous order-level rules)", file=sys.stderr)
    print(f"[guess_host] host_final distribution:\n"
          f"{df['host_final'].value_counts().to_string()}", file=sys.stderr)
    print(f"[guess_host] written: {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
