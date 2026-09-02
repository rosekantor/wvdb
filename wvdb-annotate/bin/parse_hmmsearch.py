#!/usr/bin/env python3
"""
parse_hmmsearch.py — parse hmmsearch --domtblout output to a clean TSV.

The domain table format has whitespace-delimited columns with a complex
header. This script parses it into a standard TSV with one row per
domain hit, filtered by i-Evalue and c-Evalue.

Output columns:
  target_name      protein ID (query sequence)
  target_accession protein accession (or -)
  target_len       protein length
  query_name       HMM profile name
  query_accession  HMM profile accession
  query_len        HMM profile length
  full_evalue      E-value for full sequence
  full_score       bit score for full sequence
  full_bias        bias for full sequence
  domain_n         domain number
  domain_of        total domains in sequence
  c_evalue         conditional E-value for this domain
  i_evalue         independent E-value for this domain
  domain_score     bit score for this domain
  domain_bias      bias for this domain
  hmm_from         HMM start position
  hmm_to           HMM end position
  ali_from         alignment start in target
  ali_to           alignment end in target
  env_from         envelope start in target
  env_to           envelope end in target
  acc              mean posterior probability
  description      target description

Usage:
  parse_hmmsearch.py \\
      --domtbl profile.hmmsearch.domtbl \\
      --out    profile.hmmsearch.tsv \\
      [--max-ievalue 1e-5]
"""

import argparse
import sys
from pathlib import Path


DOMTBL_COLS = [
    'target_name', 'target_accession', 'target_len',
    'query_name', 'query_accession', 'query_len',
    'full_evalue', 'full_score', 'full_bias',
    'domain_n', 'domain_of', 'c_evalue', 'i_evalue',
    'domain_score', 'domain_bias',
    'hmm_from', 'hmm_to',
    'ali_from', 'ali_to',
    'env_from', 'env_to',
    'acc', 'description'
]

# Number of whitespace-delimited fixed fields before the description
N_FIXED = 22


def parse_domtbl(domtbl_path, max_ievalue=None):
    """Parse hmmsearch --domtblout file into a list of dicts."""
    rows = []
    with open(domtbl_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < N_FIXED:
                print(f"Warning: skipping short line: {repr(line[:80])}",
                      file=sys.stderr)
                continue
            fixed = parts[:N_FIXED]
            description = ' '.join(parts[N_FIXED:]) if len(parts) > N_FIXED else ''
            row = dict(zip(DOMTBL_COLS[:N_FIXED], fixed))
            row['description'] = description

            # Apply i-Evalue filter if requested
            if max_ievalue is not None:
                try:
                    if float(row['i_evalue']) > max_ievalue:
                        continue
                except ValueError:
                    pass

            rows.append(row)
    return rows


def main():
    p = argparse.ArgumentParser(
        description="Parse hmmsearch --domtblout to a clean TSV."
    )
    p.add_argument('--domtbl',       required=True, help='hmmsearch --domtblout file')
    p.add_argument('--out',          required=True, help='Output TSV path')
    p.add_argument('--max-ievalue',  type=float, default=None,
                   help='Maximum independent E-value to keep (default: no filter, '
                        'rely on hmmsearch --cut_tc or --domE)')
    args = p.parse_args()

    if not Path(args.domtbl).exists():
        print(f"ERROR: domtbl file not found: {args.domtbl}", file=sys.stderr)
        sys.exit(1)

    rows = parse_domtbl(args.domtbl, args.max_ievalue)

    if not rows:
        print(f"Warning: no domain hits found in {args.domtbl}", file=sys.stderr)
        # Write empty TSV with header so downstream tools don't fail
        with open(args.out, 'w') as f:
            f.write('\t'.join(DOMTBL_COLS) + '\n')
        return

    with open(args.out, 'w') as f:
        f.write('\t'.join(DOMTBL_COLS) + '\n')
        for row in rows:
            f.write('\t'.join(row.get(c, '') for c in DOMTBL_COLS) + '\n')

    print(f"parse_hmmsearch: {len(rows)} domain hits → {args.out}", file=sys.stderr)


if __name__ == '__main__':
    main()
