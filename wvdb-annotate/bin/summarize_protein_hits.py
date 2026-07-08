#!/usr/bin/env python3
"""
summarize_protein_hits.py — summarize protein search results per contig.

Reads predicted proteins from the geNomad .faa file to get total protein
counts per contig, then counts proteins with at least one hit per database
(DIAMOND or hmmsearch). Writes a per-contig summary TSV for incorporation
into merged_annotations.tsv.

Protein IDs in the .faa file are expected to follow Prodigal/geNomad naming:
  <contig_id>_<protein_number>   e.g. NODE_1_length_50000_cov_10_1

Output columns:
  contig                  contig ID
  n_proteins_total        total predicted proteins on this contig
  n_proteins_<db>_hit     proteins with ≥1 hit in each database (one col per db)

Usage:
  summarize_protein_hits.py \\
      --faa         proteins.faa \\
      --diamond-dir .              (directory containing *.diamond.tsv files) \\
      --hmmsearch-dir .            (directory containing *.hmmsearch.tsv files) \\
      --out         protein_summary.tsv
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


def protein_to_contig(protein_id):
    """
    Extract contig ID from a Prodigal/geNomad protein ID.
    Convention: <contig_id>_<integer> where the trailing integer is the
    ORF number. Handles contigs that themselves contain underscores.
    """
    # Strip trailing _N where N is an integer
    match = re.match(r'^(.+)_(\d+)$', protein_id)
    if match:
        return match.group(1)
    return protein_id


def count_proteins_per_contig(faa_path):
    """Return {contig: n_proteins} from a FASTA of predicted proteins."""
    counts = defaultdict(int)
    with open(faa_path) as f:
        for line in f:
            if line.startswith('>'):
                protein_id = line[1:].split()[0]
                contig = protein_to_contig(protein_id)
                counts[contig] += 1
    return dict(counts)


def count_hits_per_contig(tsv_path, db_name, qseqid_col=0):
    """
    Count proteins with ≥1 hit per contig from a TSV file.
    qseqid_col: column index of the query sequence ID (default 0).
    """
    p = Path(tsv_path)
    if not p.exists() or p.stat().st_size == 0:
        return {}

    hits = defaultdict(set)
    with open(p) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) <= qseqid_col:
                continue
            protein_id = parts[qseqid_col].strip()
            if protein_id in ('target_name', 'qseqid'):
                continue  # skip header if present
            contig = protein_to_contig(protein_id)
            hits[contig].add(protein_id)

    return {contig: len(proteins) for contig, proteins in hits.items()}


def main():
    p = argparse.ArgumentParser(
        description="Summarize protein search results per contig."
    )
    p.add_argument('--faa',           required=True,
                   help='geNomad predicted proteins .faa file')
    p.add_argument('--diamond-dir',   default=None,
                   help='Directory containing *.diamond.tsv files')
    p.add_argument('--hmmsearch-dir', default=None,
                   help='Directory containing *.hmmsearch.tsv files')
    p.add_argument('--out',           required=True,
                   help='Output per-contig protein summary TSV')
    args = p.parse_args()

    # Count total proteins per contig
    print(f"Counting proteins from: {args.faa}", file=sys.stderr)
    total_counts = count_proteins_per_contig(args.faa)
    df = pd.DataFrame([
        {'contig': c, 'n_proteins_total': n}
        for c, n in total_counts.items()
    ])
    print(f"  {len(df)} contigs with predicted proteins", file=sys.stderr)

    # Count hits per contig for each DIAMOND database
    if args.diamond_dir:
        diamond_dir = Path(args.diamond_dir)
        diamond_tsvs = sorted(diamond_dir.rglob('*.diamond.tsv'))
        for tsv in diamond_tsvs:
            db_name = tsv.name.replace('.diamond.tsv', '')
            col = f'n_proteins_{db_name}_hit'
            print(f"  DIAMOND hits: {db_name} ({tsv.name})", file=sys.stderr)
            # DIAMOND outfmt 6: qseqid is col 0
            hits = count_hits_per_contig(tsv, db_name, qseqid_col=0)
            df[col] = df['contig'].map(hits).fillna(0).astype(int)

    # Count hits per contig for each hmmsearch profile
    if args.hmmsearch_dir:
        hmm_dir = Path(args.hmmsearch_dir)
        hmm_tsvs = sorted(hmm_dir.rglob('*.hmmsearch.tsv'))
        for tsv in hmm_tsvs:
            profile_name = tsv.name.replace('.hmmsearch.tsv', '')
            col = f'n_proteins_{profile_name}_hit'
            print(f"  hmmsearch hits: {profile_name} ({tsv.name})", file=sys.stderr)
            # parse_hmmsearch.py TSV: target_name is col 0
            hits = count_hits_per_contig(tsv, profile_name, qseqid_col=0)
            df[col] = df['contig'].map(hits).fillna(0).astype(int)

    df.to_csv(args.out, sep='\t', index=False)
    print(f"protein_summary: {len(df)} contigs → {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
