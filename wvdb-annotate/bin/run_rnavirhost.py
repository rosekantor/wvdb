#!/usr/bin/env python3
"""
run_rnavirhost.py — generate consensus viral order classifications and run RNAVirHost.

Builds the order classification CSV required by RNAVirHost from RdRPCATCH and
geNomad outputs (RdRPCATCH order preferred; geNomad fallback; 'Unclassified' if
neither). Runs `rnavirhost predict` and reports the result CSV path.

The input FASTA is the canonical source of which contigs exist — RNAVirHost
requires orders.csv to have exactly one row per FASTA sequence, in the same
set. RdRPCATCH and geNomad outputs are left-joined onto that canonical list
(never outer-joined), so any stray or duplicate contig IDs in their output
can never inflate or shrink the row count relative to the FASTA.

Usage:
  run_rnavirhost.py \\
      -g genomad_out/<prefix>_summary/<prefix>_virus_summary.tsv \\
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

import pandas as pd
from Bio import SeqIO

from genomad_utils import load_genomad_virus_summary


# ---------------------------------------------------------------------------
# Loaders (parallel to merge_annotations.py — kept lightweight here)
# ---------------------------------------------------------------------------

def load_genomad(genomad):
    df = load_genomad_virus_summary(genomad)
    return df[['contig', 'Order']]


def load_rdrpcatch(rdrpcatch_file):
    df = pd.read_csv(rdrpcatch_file, sep='\t')
    df = df.rename(columns={'Contig_name': 'contig'})
    tax = df['MMseqs_Taxonomy_2bLCA'].fillna('').str.replace(' ', '_', regex=False)
    df['Order'] = tax.str.extract(r';o_([^;]+);?', expand=False)
    return df[['contig', 'Order']].copy()


def load_fasta_ids(fasta_path):
    """Return the list of sequence IDs in a FASTA, in file order."""
    return [rec.id for rec in SeqIO.parse(fasta_path, 'fasta')]


# ---------------------------------------------------------------------------
# Consensus order + orders.csv generation
# ---------------------------------------------------------------------------

def get_consensus_order(fasta_ids, rdrpcatch_df, genomad_df, orders_file):
    """
    Build orders.csv for RNAVirHost from RdRPCATCH + geNomad order calls.

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

    # Canonical contig list from the FASTA — left-join everything onto this
    orders_df = pd.DataFrame({'contig': fasta_ids})
    orders_df = orders_df.merge(rdrp_orders,    on='contig', how='left')
    orders_df = orders_df.merge(genomad_orders, on='contig', how='left')

    orders_df['Order_consensus'] = orders_df['Order_R']
    orders_df.loc[orders_df['Order_consensus'].isna(), 'Order_consensus'] = \
        orders_df['Order_G']
    orders_df['Order_consensus'] = orders_df['Order_consensus'].fillna('Unclassified')

    # Sanity check: row count must exactly equal FASTA sequence count
    assert len(orders_df) == len(fasta_ids), (
        f"orders_df has {len(orders_df)} rows but FASTA has {len(fasta_ids)} "
        f"sequences — this should never happen since orders_df is built "
        f"directly from fasta_ids"
    )

    # Write in the format RNAVirHost expects: unnamed first col, 'y|virus order' second col
    out_df = orders_df[['contig', 'Order_consensus']].rename(
        columns={'contig': '', 'Order_consensus': 'y|virus order'}
    )
    out_df.to_csv(orders_file, index=False)

    n_classified = (orders_df['Order_consensus'] != 'Unclassified').sum()
    n_rdrp_unmatched    = (~rdrp_orders['contig'].isin(fasta_ids)).sum()
    n_genomad_unmatched = (~genomad_orders['contig'].isin(fasta_ids)).sum()

    print(
        f"[run_rnavirhost] orders file: {len(orders_df)} contigs "
        f"(matches FASTA exactly), {n_classified} classified, "
        f"{len(orders_df) - n_classified} Unclassified",
        file=sys.stderr
    )
    if n_rdrp_unmatched > 0:
        print(f"[run_rnavirhost] warning: {n_rdrp_unmatched} RdRPCATCH contig "
              f"IDs not found in FASTA — ignored", file=sys.stderr)
    if n_genomad_unmatched > 0:
        print(f"[run_rnavirhost] warning: {n_genomad_unmatched} geNomad contig "
              f"IDs not found in FASTA — ignored", file=sys.stderr)


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
    parser.add_argument('-r', '--rdrpcatch', required=True,
                        help='RdRPCATCH annotated output TSV (same as --rdrpcatch in merge_annotations.py)')
    parser.add_argument('-f', '--fasta', required=True,
                        help='Input vOTU FASTA — canonical source of contig IDs')
    parser.add_argument('-O', '--orders-file', required=True,
                        help='Output consensus orders CSV for RNAVirHost')
    parser.add_argument('-o', '--rnavirhost-outdir', required=True,
                        help='Output directory for RNAVirHost')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing rnavirhost output directory')

    args = parser.parse_args()

    fasta_ids  = load_fasta_ids(args.fasta)
    genomad_df = load_genomad(args.genomad)
    rdrp_df    = load_rdrpcatch(args.rdrpcatch)

    get_consensus_order(fasta_ids, rdrp_df, genomad_df, args.orders_file)

    run_rnavirhost(
        fasta       = args.fasta,
        orders_file = args.orders_file,
        outdir      = args.rnavirhost_outdir,
        force       = args.force,
    )


if __name__ == '__main__':
    main()
