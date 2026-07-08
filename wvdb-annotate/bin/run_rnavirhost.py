#!/usr/bin/env python3
"""
run_rnavirhost.py — generate consensus viral order classifications and run RNAVirHost.

Builds the order classification CSV required by RNAVirHost from RdRPCATCH and
geNomad outputs (RdRPCATCH order preferred; geNomad fallback; 'Unclassified' if
neither). Runs `rnavirhost predict` and reports the result CSV path.

Usage:
  run_rnavirhost.py \\
      -g genomad_out/<prefix>_summary/<prefix>_virus_summary.tsv \\
      -c checkv_out/quality_summary.tsv \\
      -r rdrpcatch_out/rdrpcatch_output_annotated.tsv \\
      -f votus.fasta \\
      -O orders.csv \\
      -o rnavirhost_out
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Loaders (parallel to merge_annotations.py — kept lightweight here)
# ---------------------------------------------------------------------------

def load_checkv(checkv):
    df = pd.read_csv(checkv, sep='\t')
    df = df.rename(columns={'contig_id': 'contig'})
    return df


def load_genomad(genomad):
    df = pd.read_csv(genomad, sep='\t')
    df = df.rename(columns={'seq_name': 'contig'})
    df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']] = (
        df['taxonomy'].str.split(';', expand=True)
    )
    df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']] = (
        df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']]
        .replace('', np.nan)
    )
    return df[['contig', 'Order']]


def load_rdrpcatch(rdrpcatch_file):
    df = pd.read_csv(rdrpcatch_file, sep='\t')
    df = df.rename(columns={'Contig_name': 'contig'})
    tax = df['MMseqs_Taxonomy_2bLCA'].fillna('').str.replace(' ', '_', regex=False)
    df['Order'] = tax.str.extract(r';o_([^;]+);?', expand=False)
    return df[['contig', 'Order']].copy()


# ---------------------------------------------------------------------------
# Consensus order + orders.csv generation
# ---------------------------------------------------------------------------

def get_consensus_order(rdrpcatch_df, genomad_df, checkv_df, orders_file):
    """
    Build orders.csv for RNAVirHost from RdRPCATCH + geNomad + CheckV contigs.

    Priority:
      1. RdRPCATCH order
      2. geNomad order
      3. 'Unclassified'

    The output CSV format required by rnavirhost predict:
      (unnamed index col)  y|virus order
      contig_id_1          Ortervirales
      ...
    """
    orders_file = Path(orders_file)
    orders_file.parent.mkdir(parents=True, exist_ok=True)

    rdrp_orders    = rdrpcatch_df.drop_duplicates('contig').rename(columns={'Order': 'Order_R'})
    genomad_orders = genomad_df.drop_duplicates('contig').rename(columns={'Order': 'Order_G'})
    checkv_contigs = checkv_df[['contig']].drop_duplicates('contig')

    orders_df = pd.merge(rdrp_orders,    genomad_orders, on='contig', how='outer')
    orders_df = pd.merge(orders_df,      checkv_contigs, on='contig', how='outer')

    orders_df['Order_consensus'] = orders_df['Order_R']
    orders_df.loc[orders_df['Order_consensus'].isna(), 'Order_consensus'] = \
        orders_df['Order_G']
    orders_df['Order_consensus'] = orders_df['Order_consensus'].fillna('Unclassified')

    # Write in the format RNAVirHost expects: unnamed first col, 'y|virus order' second col
    out_df = orders_df[['contig', 'Order_consensus']].rename(
        columns={'contig': '', 'Order_consensus': 'y|virus order'}
    )
    out_df.to_csv(orders_file, index=False)

    n_classified = (orders_df['Order_consensus'] != 'Unclassified').sum()
    print(
        f"[run_rnavirhost] orders file: {len(orders_df)} contigs, "
        f"{n_classified} classified, "
        f"{len(orders_df) - n_classified} Unclassified",
        file=sys.stderr
    )

    return orders_df[['contig', 'Order_R', 'Order_G', 'Order_consensus']]


# ---------------------------------------------------------------------------
# RNAVirHost runner
# ---------------------------------------------------------------------------

def run_rnavirhost(fasta, orders_file, outdir, force=False):
    """
    Run `rnavirhost predict` and return the result CSV path.

    Parameters
    ----------
    fasta : str | Path
    orders_file : str | Path
    outdir : str | Path
    force : bool
        If True, remove existing outdir before running. Required for
        Nextflow -resume compatibility since work dirs may be reused.
    """
    fasta       = Path(fasta)
    orders_file = Path(orders_file)
    outdir      = Path(outdir)

    if outdir.exists():
        if force:
            print(f"  Removing existing output directory: {outdir}", file=sys.stderr)
            shutil.rmtree(outdir)
        else:
            raise FileExistsError(
                f"Output directory already exists: {outdir}\n"
                f"Use --force to overwrite."
            )

    cmd = [
        "rnavirhost", "predict",
        "-i",      str(fasta),
        "--taxa",  str(orders_file),
        "-o",      str(outdir),
    ]
    print(f"[run_rnavirhost] running: {' '.join(cmd)}", file=sys.stderr)

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout, file=sys.stderr)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    if result.returncode != 0:
        raise RuntimeError(
            f"rnavirhost predict failed (exit {result.returncode})"
        )

    # Try both output path conventions across rnavirhost versions
    csv_path = outdir / "predict" / "result.csv"
    if not csv_path.exists():
        csv_path = outdir / "result.csv"
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Expected RNAVirHost output not found: {csv_path}\n"
            f"Check rnavirhost stdout/stderr above."
        )

    print(f"[run_rnavirhost] result: {csv_path}", file=sys.stderr)
    return csv_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate consensus viral order classifications from geNomad and "
            "RdRPCATCH, then run RNAVirHost host prediction."
        )
    )
    parser.add_argument('-g', '--genomad', required=True,
                        help='geNomad <prefix>_virus_summary.tsv (same as --genomad in merge_annotations.py)')
    parser.add_argument('-c', '--checkv', required=True,
                        help='CheckV quality_summary.tsv (same as --checkv in merge_annotations.py)')
    parser.add_argument('-r', '--rdrpcatch', required=True,
                        help='RdRPCATCH annotated output TSV (same as --rdrpcatch in merge_annotations.py)')
    parser.add_argument('-f', '--fasta', required=True,
                        help='Input vOTU FASTA')
    parser.add_argument('-O', '--orders-file', required=True,
                        help='Output consensus orders CSV for RNAVirHost')
    parser.add_argument('-o', '--rnavirhost-outdir', required=True,
                        help='Output directory for RNAVirHost')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing rnavirhost output directory')

    args = parser.parse_args()

    checkv_df  = load_checkv(args.checkv)
    genomad_df = load_genomad(args.genomad)
    rdrp_df    = load_rdrpcatch(args.rdrpcatch)

    get_consensus_order(rdrp_df, genomad_df, checkv_df, args.orders_file)

    run_rnavirhost(
        fasta       = args.fasta,
        orders_file = args.orders_file,
        outdir      = args.rnavirhost_outdir,
        force       = args.force,
    )


if __name__ == '__main__':
    main()
