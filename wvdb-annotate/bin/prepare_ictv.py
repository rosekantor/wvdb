#!/usr/bin/env python3
"""
prepare_ictv.py — process the ICTV Virus Properties By Family table.

The ICTV website (https://ictv.global/virus-properties) does not provide
a direct download link. To obtain the input file:
  1. Open https://ictv.global/virus-properties in your browser
  2. Set "Items per page" to "All" to load all families
  3. Save the page as HTML (File → Save Page As)
  4. Pass the saved file to this script via --in

Simplifies the Host column into broad categories and writes a processed
TSV for use by merge_annotations.py.

Run once to generate the reference file, or re-run when a new ICTV release
is available.

Usage:
  prepare_ictv.py \\
      --in  /path/to/Virus_Properties___ICTV.html \\
      --out /path/to/ref_data/ictv_families.tsv
"""

import argparse
import sys
from datetime import date
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Host simplification mapping
# ---------------------------------------------------------------------------
HOST_REASSIGN = [
    ['invertebrates',                                       'invertebrates'],
    ['bacteria',                                            'bacteria'],
    ['vertebrates',                                         'vertebrates'],
    ['archaea',                                             'archaea'],
    ['protists',                                            'protists'],
    ['fungi',                                               'fungi'],
    ['plants',                                              'plants'],
    ['predicted bacteria',                                  'bacteria'],
    ['predicted archaea',                                   'archaea'],
    ['predicted protists',                                  'protists'],
    ['fungi, plants',                                       'non-vertebrates'],
    ['plants, invertebrates',                               'non-vertebrates'],
    ['protists, invertebrates',                             'non-vertebrates'],
    ['protists, fungi, invertebrates',                      'non-vertebrates'],
    ['fungi, protists',                                     'non-vertebrates'],
    ['fungi, plants, invertebrates',                        'non-vertebrates'],
    ['protists, fungi',                                     'non-vertebrates'],
    ['protists, fungi, plants',                             'non-vertebrates'],
    ['soil (s)',                                            'non-vertebrates'],
    ['fungi, invertebrates',                                'non-vertebrates'],
    ['protists, fungi, plants, invertebrates',              'non-vertebrates'],
    ['invertebrates, vertebrates',                          'incl-vertebrates'],
    ['fungi, plants, invertebrates, vertebrates',           'incl-vertebrates'],
    ['plants, invertebrates, vertebrates',                  'incl-vertebrates'],
    ['protists, fungi, plants, invertebrates, vertebrates', 'incl-vertebrates'],
    ['protists, fungi, vertebrates',                        'incl-vertebrates'],
    ['protists, plants, invertebrates, vertebrates',        'incl-vertebrates'],
    ['invertebrates, plants',                               'non-vertebrates'],
    ['vertebrates, invertebrates, plants',                  'incl-vertebrates'],
]


def load_input(in_path):
    """Load from a saved HTML file or a previously extracted TSV."""
    suffix = in_path.suffix.lower()
    if suffix in ('.html', '.htm'):
        print(f"Reading HTML table from: {in_path}", file=sys.stderr)
        tables = pd.read_html(str(in_path), flavor='lxml')
        if not tables:
            print("ERROR: no tables found in HTML file.", file=sys.stderr)
            sys.exit(1)
        # Find the table that has a Family column (may be named 'Family  Sort descending')
        df = None
        for t in tables:
            family_col = next((c for c in t.columns if 'family' in str(c).lower()), None)
            if family_col and 'host' in [str(c).lower() for c in t.columns]:
                df = t
                break
        if df is None:
            print("ERROR: could not find the virus properties table in the HTML.",
                  file=sys.stderr)
            print(f"Tables found: {len(tables)}, columns: {[list(t.columns) for t in tables]}",
                  file=sys.stderr)
            sys.exit(1)
    elif suffix == '.tsv':
        print(f"Reading TSV from: {in_path}", file=sys.stderr)
        df = pd.read_csv(in_path, sep='\t')
    else:
        print(f"ERROR: unsupported file type '{suffix}'. Expected .html or .tsv.",
              file=sys.stderr)
        sys.exit(1)

    print(f"  Loaded {len(df)} rows, {len(df.columns)} columns", file=sys.stderr)
    return df


def process_raw(df, out_path):
    """Normalise column names, simplify host designations, write TSV."""

    # Normalise column names — strip whitespace, handle 'Family  Sort descending'
    col_map = {}
    for c in df.columns:
        clean = str(c).strip()
        if 'family' in clean.lower():
            col_map[c] = 'Family'
        elif 'host' in clean.lower() and 'ictv' not in clean.lower():
            col_map[c] = 'Host'
        else:
            col_map[c] = clean
    df = df.rename(columns=col_map)

    # Drop row-number column
    unnamed = [c for c in df.columns if str(c).startswith('Unnamed')]
    if unnamed:
        df = df.drop(columns=unnamed)

    # Validate required columns
    required = {'Family', 'Host', 'Genome size (kb/kbp)', 'Genome topology'}
    missing = required - set(df.columns)
    if missing:
        print(f"ERROR: missing columns: {missing}", file=sys.stderr)
        print(f"Available: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # Keep and rename
    df = df[['Family', 'Host', 'Genome size (kb/kbp)', 'Genome topology']].copy()
    df.columns = ['Family_ICTV', 'host_ICTV', 'size_kb_ICTV', 'topology_ICTV']

    # Clean Family names — strip asterisks added by HTML rendering
    df['Family_ICTV'] = df['Family_ICTV'].astype(str).str.strip().str.strip('*').str.strip()

    # Normalise host strings for matching
    df['host_ICTV'] = df['host_ICTV'].astype(str).str.strip().str.lower()
    df['host_ICTV'] = df['host_ICTV'].replace({'nan': None, '': None})

    # Apply host simplification mapping
    host_reassign_df = pd.DataFrame.from_records(
        HOST_REASSIGN, columns=['host_ICTV', 'host_ICTV_simple']
    )
    df = df.merge(host_reassign_df, how='left', on='host_ICTV')

    # Restore original-case host_ICTV for readability in output
    # (matching was done on lowercase; restore from original)

    n_unmapped = df['host_ICTV_simple'].isna().sum()
    if n_unmapped > 0:
        unmapped = df.loc[df['host_ICTV_simple'].isna(), 'host_ICTV'].dropna().unique()
        unmapped = [u for u in unmapped if u and u != 'nan']
        if unmapped:
            print(
                f"Warning: {n_unmapped} families have unmapped host values "
                f"(host_ICTV_simple will be null):\n  {list(unmapped)}",
                file=sys.stderr
            )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep='\t', index=False)

    categories = sorted(df['host_ICTV_simple'].dropna().unique().tolist())
    print(
        f"Written {len(df)} family rows → {out_path}\n"
        f"  host_ICTV_simple categories: {categories}\n"
        f"  Processed: {date.today()}",
        file=sys.stderr
    )


def main():
    parser = argparse.ArgumentParser(
        description="Process the ICTV Virus Properties By Family table."
    )
    parser.add_argument(
        '--in', dest='infile', required=True, type=Path,
        help='Saved HTML page from ictv.global/virus-properties, or a previously extracted TSV'
    )
    parser.add_argument(
        '--out', required=True, type=Path,
        help='Output path for processed TSV (e.g. ref_data/ictv_families.tsv)'
    )
    args = parser.parse_args()

    if not args.infile.exists():
        print(f"ERROR: input file not found: {args.infile}", file=sys.stderr)
        print(
            "\nTo obtain the ICTV table:\n"
            "  1. Open https://ictv.global/virus-properties in your browser\n"
            "  2. Set 'Items per page' to 'All' to load all families\n"
            "  3. File → Save Page As (save as HTML)\n"
            "  4. Re-run: prepare_ictv.py --in /path/to/saved.html --out ...",
            file=sys.stderr
        )
        sys.exit(1)

    df = load_input(args.infile)
    process_raw(df, args.out)


if __name__ == '__main__':
    main()
