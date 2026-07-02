import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess

def load_checkv(checkv):
    checkv_df = pd.read_csv(checkv, sep='\t')
    checkv_df = checkv_df.rename(columns={'contig_id': 'contig'})
    return checkv_df


def load_genomad(genomad):
    genomad_tax_df = pd.read_csv(genomad, sep='\t')
    genomad_tax_df = genomad_tax_df.rename(columns={'seq_name': 'contig'})
    genomad_tax_df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']] = (
        genomad_tax_df['taxonomy'].str.split(';', expand=True)
    )
    genomad_tax_df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']] = (
        genomad_tax_df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']].replace('', np.nan)
    )
    return genomad_tax_df[['contig', 'Order']]


def load_rdrpcatch(rdrpcatch_file):
    rdrp_df = pd.read_csv(rdrpcatch_file, sep='\t')
    rdrp_df = rdrp_df.rename(columns={'Contig_name': 'contig'})
    
    tax = rdrp_df['MMseqs_Taxonomy_2bLCA'].fillna('').str.replace(' ', '_', regex=False)
    rdrp_df['Order'] = tax.str.extract(r';o_([^;]+);?', expand=False)

    rdrp_tax_df = rdrp_df[['contig', 'Order']].copy()

    return rdrp_tax_df


def get_consensus_order(rdrpcatch_df, genomad_df, checkv_df, orders_file):
    """
    Make orders.csv file for RNAVirHost.

    Priority:
    1. RdRpCATCH order
    2. geNomad order
    3. Unclassified
    """
    orders_file = Path(orders_file)
    orders_file.parent.mkdir(parents=True, exist_ok=True)

    rdrp_orders = rdrpcatch_df.drop_duplicates('contig').rename(columns={'Order': 'Order_R'})
    genomad_orders = genomad_df.drop_duplicates('contig').rename(columns={'Order': 'Order_G'})
    # make sure all contigs exist in the dataframe (they all exist in checkv but not all necessarily exist in genomad or rdrpcatch)
    checkv_contigs = checkv_df[['contig']].drop_duplicates('contig')

    orders_df = pd.merge(rdrp_orders, genomad_orders, on='contig', how='outer')
    orders_df = pd.merge(orders_df, checkv_contigs, on='contig', how='outer')

    orders_df['Order_consensus'] = orders_df['Order_R']
    orders_df.loc[orders_df['Order_consensus'].isna(), 'Order_consensus'] = orders_df['Order_G']
    orders_df['Order_consensus'] = orders_df['Order_consensus'].fillna('Unclassified')
    out_df = orders_df[['contig', 'Order_consensus']].rename(
        columns={'contig': '', 'Order_consensus': 'y|virus order'}
    )

    out_df.to_csv(orders_file, index=False)

    return orders_df[['contig', 'Order_R', 'Order_G', 'Order_consensus']]


def run_rnavirhost(fasta, orders_file, outdir):
    """
    Run rnavirhost predict and return the expected result CSV plus exit code.
    """
    fasta = Path(fasta)
    orders_file = Path(orders_file)
    outdir = Path(outdir)

    if outdir.exists():
        raise FileExistsError(f"Output directory already exists: {outdir}")

    cmd = [
        "rnavirhost",
        "predict",
        "-i", str(fasta),
        "--taxa", str(orders_file),
        "-o", str(outdir),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    csv_path = outdir / "predict" / "result.csv"

    if result.returncode != 0:
        raise RuntimeError(
            f"rnavirhost failed with exit code {result.returncode}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    if not csv_path.exists():
        raise FileNotFoundError(f"Expected output was not found: {csv_path}")

    return csv_path, result.returncode


def main():
    parser = argparse.ArgumentParser(
    description=(
    "Use RNAVirHost with custom Order-level classifications." \
    "Get consensus orders from geNomad and RdRpCATCH, feed to RNAVirHost, run RNAVirHost.")
    )
    parser.add_argument(
        '-g', '--genomad_virus_summary', required=True, type=str,
        help='geNomad virus_summary.tsv'
    )
    parser.add_argument(
        '-c', '--checkv_tsv', required=True, type=str,
        help='CheckV quality_summary.tsv'
    )
    parser.add_argument(
        '-r', '--rdrpcatch_file', required=True, type=str,
        help='RdRpCATCH taxonomy tsv'
    )
    parser.add_argument(
        '-f', '--fasta', required=True, type=str,
        help='fasta file of viral genomes'
    )
    parser.add_argument(
        '-O', '--orders_file', required=True, type=str,
        help='Output consensus order file'
    )
    parser.add_argument(
        '-o', '--rnavirhost_outdir', required=True, type=str,
        help='Output directory for RNAVirHost'
    )

    args = parser.parse_args()

    checkv_df = load_checkv(args.checkv_tsv)
    genomad_tax_df = load_genomad(args.genomad_summary)
    rdrp_tax_df = load_rdrpcatch(args.rdrpcatch_file)

    # combine files
    get_consensus_order(rdrp_tax_df, genomad_tax_df, checkv_df, args.orders_file)
    rnavirhost_csv, exitcode = run_rnavirhost(args.fasta, args.orders_file, args.rnavirhost_outdir)

    rnavirhost_csv = Path(args.rnavirhost_outdir) / "predict" / "result.csv"
    if rnavirhost_csv.exists() and rnavirhost_csv.stat().st_size > 0:
        print(rnavirhost_csv)
    else:
        rnavirhost_csv, _ = run_rnavirhost(args.fasta, args.orders_file, args.rnavirhost_outdir)


if __name__ == '__main__':
    main()
