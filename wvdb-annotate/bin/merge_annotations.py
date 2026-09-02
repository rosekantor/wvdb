#!/usr/bin/env python3
"""
merge_annotations.py — build a merged per-vOTU annotation table for WVDB.

Combines outputs from CheckV, geNomad, RdRPCATCH (optional), BLASTn against
multiple databases, Entrez host/lineage lookup (core_nt only), RNAVirHost
(optional), and the ICTV family reference table.

BLASTn results are read from a directory containing one ANI TSV per database
(<db_name>.blastn.ani.tsv), produced by blastani_nayfach.py. All databases
contribute columns to the merged table; core_nt additionally triggers an
Entrez lookup for host and lineage information.

Usage:
  merge_annotations.py \\
      --checkv        checkv_out/quality_summary.tsv \\
      --genomad       genomad_out/<prefix>_summary/<prefix>_virus_summary.tsv \\
      --ictv-fam      ref_data/ictv_families.tsv \\
      --entrez-email  user@institution.edu \\
      --out           merged_annotations.tsv \\
      [--blastn-dir   blastn/] \\
      [--rdrpcatch    rdrpcatch_out/rdrpcatch_output_annotated.tsv] \\
      [--rnavirhost   rnavirhost_out/predict/result.csv] \\
      [--host-lineage-cache  host_lineage.tsv]
"""

import argparse
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO

from genomad_utils import load_genomad_virus_summary


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_checkv(path):
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={'contig_id': 'contig'})
    return df


def load_genomad(path):
    return load_genomad_virus_summary(path)


def load_rdrpcatch(path):
    """Load RdRPCATCH output and parse MMseqs taxonomy string."""
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={'Contig_name': 'contig'})
    df['MMseqs_Taxonomy_2bLCA'] = (
        df['MMseqs_Taxonomy_2bLCA'].fillna('').str.replace(' ', '_', regex=False)
    )
    tax = df['MMseqs_Taxonomy_2bLCA']
    df['Domain']  = tax.str.extract(r'd_([^;]+);?')
    df['Realm']   = tax.str.extract(r'd_[^;]+;-_([^;]+);?')
    df['Kingdom'] = tax.str.extract(r';k_([^;]+);?')
    df['Phylum']  = tax.str.extract(r';p_([^;]+);?')
    df['Class']   = tax.str.extract(r';c_([^;]+);?')
    df['Order']   = tax.str.extract(r';o_([^;]+);?')
    df['Family']  = tax.str.extract(r';f_([^;]+);?')
    df['Genus']   = tax.str.extract(r';g_([^;]+);?')
    df['Species'] = tax.str.extract(r';s_([^;]+);?')
    df['Lowest']  = tax.str.extract(r'_([^;]+)$')

    return df[[
        'contig', 'Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order',
        'Family', 'Genus', 'Species', 'Lowest', 'MMseqs_Taxonomy_2bLCA',
        'Best_hit_bitscore'
    ]].copy()


def load_rnavirhost(path):
    df = pd.read_csv(
        path, skiprows=1,
        names=['contig', 'Order', 'host_domain', 'host_chordata', 'evidence']
    )
    df = df.replace(np.nan, 'unknown')
    df = df.rename(columns={
        'host_domain':    'rnavirhost_domain',
        'host_chordata':  'rnavirhost_chordata',
    })
    return df[['contig', 'rnavirhost_domain', 'rnavirhost_chordata']]


def process_blastn(tsv_path, db_name):
    """
    Load a blastn ANI TSV (output of blastani_nayfach.py), classify hits as
    genus or species level, keep best hit per query, and rename columns to
    include the database name.
    """
    df = pd.read_csv(tsv_path, sep='\t')
    df['hit_level'] = 'none'
    df.loc[(df.pid >= 70) & (df.qcov >= 85), 'hit_level'] = 'genus'
    df.loc[(df.pid >= 95) & (df.qcov >= 85), 'hit_level'] = 'species'
    df['idXcov'] = df.pid * df.qcov
    idx_max = df.groupby('qname')['idXcov'].idxmax()
    df = df.loc[idx_max]
    df = df.rename(columns={
        'qname':      'contig',
        'tname':      f'tname_{db_name}',
        'pid':        f'pid_{db_name}',
        'qcov':       f'qcov_{db_name}',
        'hit_level':  f'hit_level_{db_name}',
    })
    df = df.drop(columns=['num_alns', 'tcov', 'idXcov'], errors='ignore')
    return df


def load_all_blastn(blastn_dir):
    """
    Discover all <db_name>.blastn.ani.tsv files in blastn_dir.
    Returns a dict {db_name: DataFrame} and separately the core_nt DataFrame
    (None if core_nt not present).
    """
    blastn_dir = Path(blastn_dir)
    results = {}
    for tsv in sorted(blastn_dir.rglob('*.blastn.ani.tsv')):
        # Derive db_name from filename: <db_name>.blastn.ani.tsv
        db_name = tsv.name.replace('.blastn.ani.tsv', '')
        print(f"  Loading BLASTn results: {db_name} ({tsv})", file=sys.stderr)
        results[db_name] = process_blastn(tsv, db_name)
    return results


# ---------------------------------------------------------------------------
# Entrez host/lineage lookup (core_nt only)
# ---------------------------------------------------------------------------

def batch_extract_hosts_lineages(genbank_ids, batch_size=100, max_retries=3):
    """
    Batch Entrez efetch for host, isolation_source, and lineage.

    The Entrez.efetch() handle is read into a string and parsed via
    io.StringIO() rather than handed directly to SeqIO.parse(). Passing
    the live handle straight to SeqIO.parse() is fragile across Biopython
    versions — Bio.File.as_handle() can fail to recognize the handle type
    Entrez.efetch() returns and incorrectly try to treat it as a file
    path, raising a TypeError deep inside Biopython rather than a clear
    network/parsing error. Reading to a string first sidesteps that
    entirely and also lets us retry a failed batch instead of crashing
    the whole run on one bad response.
    """
    results = []
    for i in range(0, len(genbank_ids), batch_size):
        batch = genbank_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        print(f"  Entrez query batch {batch_num} "
              f"({len(batch)} IDs)", file=sys.stderr)

        data = None
        for attempt in range(1, max_retries + 1):
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=",".join(batch),
                    rettype="gb", retmode="text"
                )
                data = handle.read()
                handle.close()
                break
            except Exception as e:
                print(f"    attempt {attempt}/{max_retries} failed: {e}",
                      file=sys.stderr)
                if attempt == max_retries:
                    print(f"  WARNING: batch {batch_num} failed after "
                          f"{max_retries} attempts — skipping "
                          f"{len(batch)} ID(s): {batch}", file=sys.stderr)

        if data is None:
            continue

        try:
            records = list(SeqIO.parse(io.StringIO(data), "genbank"))
        except Exception as e:
            print(f"  WARNING: failed to parse batch {batch_num} response "
                  f"({e}) — skipping {len(batch)} ID(s): {batch}",
                  file=sys.stderr)
            continue

        parsed_ids = {rec.id for rec in records}
        missing = set(batch) - parsed_ids
        if missing:
            print(f"  WARNING: {len(missing)} ID(s) in batch {batch_num} "
                  f"had no GenBank record returned: {sorted(missing)}",
                  file=sys.stderr)

        for rec in records:
            hosts, isolation_source = [], None
            for feat in rec.features:
                if feat.type == "source":
                    hosts.extend(feat.qualifiers.get("host", []))
                    iso = feat.qualifiers.get("isolation_source", [None])[0]
                    if iso:
                        isolation_source = iso
            lineage = rec.annotations.get("taxonomy", [])
            results.append((
                rec.id,
                ','.join(hosts) if hosts else None,
                isolation_source,
                ';'.join(lineage) if lineage else None,
            ))
    return results


def get_host_lineage(corent_df, cache_path, target_col='tname_nt'):
    """
    Look up host and lineage for core_nt BLAST hits via Entrez.

    Uses incremental caching, not all-or-nothing: any target_col values
    already present in the cache WITH a populated hosts_ntBlastHit are
    reused; everything else (never queried, or queried previously but
    left null due to a failed/skipped batch) is retried via Entrez. This
    mirrors the caching pattern in categorize_host_terms.py.

    An earlier version of this function trusted any existing non-empty
    cache file unconditionally and returned it as-is — if a prior run's
    Entrez batches failed (contig/target_col get populated by a merge
    step regardless of Entrez success, but the data columns stay null),
    every subsequent run silently reused that broken cache forever, and
    had no way to pick up contigs from a newer pipeline run. This
    rewrite avoids both problems.
    """
    cache_cols = ['contig', target_col, 'hosts_ntBlastHit',
                  'isolation_source_ntBlastHit', 'lineage_ntBlastHit',
                  'lowest_tax_ntBlastHit']

    existing_cache = pd.DataFrame(columns=cache_cols)
    if cache_path and Path(cache_path).is_file() and Path(cache_path).stat().st_size > 0:
        print(f"  Loading Entrez cache: {cache_path}", file=sys.stderr)
        loaded = pd.read_csv(cache_path, sep='\t')
        missing_cols = set(cache_cols) - set(loaded.columns)
        if missing_cols:
            print(f"  Warning: cache missing column(s) {missing_cols} — "
                  f"treating cache as empty", file=sys.stderr)
        else:
            existing_cache = loaded

    # Derive hit_level column name from the tname column (e.g. tname_nt → hit_level_nt,
    # tname_core_nt → hit_level_core_nt) — avoids hardcoding 'nt'
    hit_level_col = target_col.replace('tname_', 'hit_level_')
    if hit_level_col not in corent_df.columns:
        # Fallback: find any hit_level column
        hit_level_col = next(
            (c for c in corent_df.columns if c.startswith('hit_level_')), None
        )
        if hit_level_col is None:
            print("  Warning: no hit_level column found — skipping Entrez lookup.",
                  file=sys.stderr)
            return pd.DataFrame(columns=cache_cols)

    for_entrez = corent_df[
        corent_df[hit_level_col].isin(['genus', 'species'])
    ].copy()

    if len(for_entrez) == 0:
        print("  No core_nt genus/species hits — skipping Entrez lookup.",
              file=sys.stderr)
        return pd.DataFrame(columns=cache_cols)

    all_ids = set(for_entrez[target_col].unique())

    # A cached ID counts as "resolved" only if it has a non-null host —
    # a row with contig/target_col populated but a null host means the
    # earlier Entrez call for that ID failed and should be retried.
    if not existing_cache.empty:
        resolved_ids = set(
            existing_cache.loc[existing_cache['hosts_ntBlastHit'].notna(), target_col]
        )
    else:
        resolved_ids = set()

    ids_to_query = sorted(all_ids - resolved_ids)

    print(f"  {len(all_ids)} unique GenBank ID(s) needed, "
          f"{len(resolved_ids & all_ids)} already resolved in cache, "
          f"{len(ids_to_query)} to query", file=sys.stderr)

    if ids_to_query:
        raw = batch_extract_hosts_lineages(ids_to_query)

        new_results = pd.DataFrame(
            raw,
            columns=[target_col, 'hosts_ntBlastHit',
                     'isolation_source_ntBlastHit', 'lineage_ntBlastHit']
        )
        new_results['lowest_tax_ntBlastHit'] = (
            new_results['lineage_ntBlastHit'].str.extract(r'.*;(.*$)')
        )

        # Any ID that was queried but got no record back still needs a
        # row (with nulls) so it isn't silently retried forever and so
        # the merge below doesn't drop it — but it also shouldn't be
        # treated as "resolved" next run, so leave hosts_ntBlastHit null.
        queried_ids_returned = set(new_results[target_col])
        missing = set(ids_to_query) - queried_ids_returned
        if missing:
            filler = pd.DataFrame({
                target_col: sorted(missing),
                'hosts_ntBlastHit': None,
                'isolation_source_ntBlastHit': None,
                'lineage_ntBlastHit': None,
                'lowest_tax_ntBlastHit': None,
            })
            new_results = pd.concat([new_results, filler], ignore_index=True)

        # Merge new results into the cache: replace any existing row for
        # the same target_col value (in case it's being retried), keep
        # everything else
        if not existing_cache.empty:
            existing_data = existing_cache.drop(columns=['contig']).drop_duplicates(target_col)
            existing_data = existing_data[~existing_data[target_col].isin(new_results[target_col])]
            updated_data = pd.concat([existing_data, new_results], ignore_index=True)
        else:
            updated_data = new_results
    else:
        updated_data = existing_cache.drop(columns=['contig']).drop_duplicates(target_col) \
            if not existing_cache.empty else pd.DataFrame(columns=cache_cols[1:])

    # Rebuild the full contig-level result for the CURRENT run's for_entrez set
    results_df = for_entrez[['contig', target_col]].merge(
        updated_data, on=target_col, how='left'
    )

    if cache_path:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(cache_path, sep='\t', index=False)
        print(f"  Entrez results cached: {cache_path}", file=sys.stderr)

    n_resolved = results_df['hosts_ntBlastHit'].notna().sum()
    print(f"  {n_resolved}/{len(results_df)} contig(s) have a resolved host "
          f"after this run", file=sys.stderr)

    return results_df


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------

def merge_master_table(
    checkv_df, genomad_df, rdrpcatch_df,
    blastn_dbs, entrez_df,
    rnavirhost_df, ictv_fam_df,
):
    """
    Merge all annotation sources into one per-vOTU table.

    Merge order:
      1. geNomad + RdRPCATCH  (outer join — taxonomy from both)
      2. + CheckV              (outer join)
      3. + RNAVirHost          (outer join, optional)
      4. + core_nt BLASTn      (left join)
      5. + Entrez host/lineage (left join on core_nt hits)
      6. + other BLASTn dbs   (left join, one per db)
      7. Consensus order/family columns
      8. + ICTV family lookup  (left join on Family_consensus)
    """
    # --- geNomad ---
    genomad_merge = genomad_df[[
        'contig', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family',
        'genomad_provirus', 'genomad_provirus_start', 'genomad_provirus_end'
    ]].drop_duplicates('contig')

    # --- RdRPCATCH ---
    if rdrpcatch_df is not None:
        rdrp_merge = rdrpcatch_df[[
            'contig', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order',
            'Family', 'Genus', 'Species', 'Lowest',
            'MMseqs_Taxonomy_2bLCA', 'Best_hit_bitscore'
        ]].drop_duplicates('contig')
        df = pd.merge(
            genomad_merge, rdrp_merge,
            on='contig', how='outer', suffixes=('_gNd', '_RdRp')
        )
        df = df.rename(columns={
            'Genus':   'Genus_RdRp',
            'Species': 'Species_RdRp',
            'Lowest':  'Lowest_RdRp',
        })
    else:
        df = genomad_merge.rename(columns={
            'Realm':   'Realm_gNd',
            'Kingdom': 'Kingdom_gNd',
            'Phylum':  'Phylum_gNd',
            'Class':   'Class_gNd',
            'Order':   'Order_gNd',
            'Family':  'Family_gNd',
        })
        for col in ['Realm_RdRp', 'Kingdom_RdRp', 'Phylum_RdRp', 'Class_RdRp',
                    'Order_RdRp', 'Family_RdRp', 'Genus_RdRp', 'Species_RdRp',
                    'Lowest_RdRp', 'MMseqs_Taxonomy_2bLCA', 'Best_hit_bitscore']:
            df[col] = np.nan

    # --- CheckV ---
    df = df.merge(checkv_df.drop_duplicates('contig'), on='contig', how='outer')

    # --- RNAVirHost ---
    if rnavirhost_df is not None:
        df = df.merge(rnavirhost_df.drop_duplicates('contig'), on='contig', how='outer')
        df['rnavirhost_domain'] = df['rnavirhost_domain'].fillna('Unknown')
    else:
        df['rnavirhost_domain'] = 'Unknown'
        df['rnavirhost_chordata'] = np.nan

    # --- core_nt BLASTn + Entrez ---
    if 'core_nt' in blastn_dbs:
        df = df.merge(blastn_dbs['core_nt'], on='contig', how='left')
        if entrez_df is not None and len(entrez_df) > 0:
            # tname column name depends on db_name: core_nt → tname_core_nt
            tname_col = next(
                (c for c in entrez_df.columns if c.startswith('tname_')), None
            )
            if tname_col and tname_col in df.columns:
                df = df.merge(entrez_df, on=['contig', tname_col], how='left')

    # --- Other BLASTn databases ---
    for db_name, blast_df in blastn_dbs.items():
        if db_name == 'core_nt':
            continue
        df = df.merge(blast_df, on='contig', how='left')

    # --- Consensus taxonomy (RdRp preferred, geNomad fallback) ---
    df['Order_consensus'] = df.get('Order_RdRp', pd.Series(dtype=str))
    df.loc[df['Order_consensus'].isna(), 'Order_consensus'] = df.get('Order_gNd')
    df['Order_consensus'] = df['Order_consensus'].fillna('Unclassified')

    df['Family_consensus'] = df.get('Family_RdRp', pd.Series(dtype=str))
    df.loc[df['Family_consensus'].isna(), 'Family_consensus'] = df.get('Family_gNd')
    df['Family_consensus'] = df['Family_consensus'].fillna('Unclassified')

    # --- ICTV family lookup ---
    df = df.merge(
        ictv_fam_df,
        how='left',
        left_on='Family_consensus',
        right_on='Family_ICTV'
    )
    df['host_ICTV_simple'] = df['host_ICTV_simple'].fillna('Unknown')

    return df


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def report(genomad_df, checkv_df, rdrpcatch_df, rnavirhost_df):
    """Print annotation count summary to stderr."""
    provirus = genomad_df[genomad_df['genomad_provirus'] == True]
    high_kmer = checkv_df[checkv_df['kmer_freq'] > 1.2]

    if len(provirus) > 0:
        print(
            f"Note: geNomad flags {len(provirus)} contig(s) as containing an "
            f"integrated provirus region that CheckV did not separately flag "
            f"as contamination — this may be a real phage/phage-plasmid with "
            f"genuine host genes (especially if reproducibly assembled across "
            f"multiple samples), not necessarily a chimeric artifact. See the "
            f"genomad_provirus column in the merged output for review: "
            f"{provirus['contig'].tolist()}",
            file=sys.stderr
        )
    if len(high_kmer) > 0:
        print(
            f"Warning: CheckV reports {len(high_kmer)} contigs with "
            f"kmer_freq > 1.2 — consider trimming: "
            f"{high_kmer['contig'].tolist()}",
            file=sys.stderr
        )

    rdrp_n = len(rdrpcatch_df) if rdrpcatch_df is not None else 'not run'
    rnavirh_n = len(rnavirhost_df) if rnavirhost_df is not None else 'not run'
    if rdrpcatch_df is not None:
        if len(rdrpcatch_df) != len(rdrpcatch_df['contig'].unique()):
            print("Warning: RdRPCATCH output has multiple RdRp genes per contig",
                  file=sys.stderr)

    print(
        f"[merge_annotations] annotation counts\n"
        f"  CheckV    : {len(checkv_df)}\n"
        f"  geNomad   : {len(genomad_df)}\n"
        f"  RdRPCATCH : {rdrp_n}\n"
        f"  RNAVirHost: {rnavirh_n}",
        file=sys.stderr
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Merge viral annotation outputs into a per-vOTU table."
    )
    # Required
    p.add_argument('--checkv',        required=True,
                   help='CheckV quality_summary.tsv')
    p.add_argument('--genomad',       required=True,
                   help='geNomad <prefix>_virus_summary.tsv')
    p.add_argument('--ictv-fam',      required=True,
                   help='Processed ICTV family TSV (ref_data/ictv_families.tsv)')
    p.add_argument('--entrez-email',  required=True,
                   help='Email for NCBI Entrez requests')
    p.add_argument('--out',           required=True,
                   help='Output merged annotation TSV')
    # Optional
    p.add_argument('--blastn-dir',    default=None,
                   help='Directory containing <db_name>.blastn.ani.tsv files')
    p.add_argument('--rdrpcatch',     default=None,
                   help='RdRPCATCH annotated output TSV')
    p.add_argument('--rnavirhost',    default=None,
                   help='RNAVirHost result.csv (rnavirhost_out/result.csv or rnavirhost_out/predict/result.csv)')
    p.add_argument('--host-lineage-cache', default=None,
                   help='Cache TSV for Entrez host/lineage results '
                        '(read if exists, write after querying)')
    return p.parse_args()


def main():
    args = parse_args()

    Entrez.email = args.entrez_email

    # --- Load required inputs ---
    print("Loading CheckV...", file=sys.stderr)
    checkv_df = load_checkv(args.checkv)

    print("Loading geNomad...", file=sys.stderr)
    genomad_df = load_genomad(args.genomad)

    print("Loading ICTV family table...", file=sys.stderr)
    ictv_fam_df = pd.read_csv(args.ictv_fam, sep='\t')

    # --- Load optional inputs ---
    rdrpcatch_df = None
    if args.rdrpcatch and not args.rdrpcatch.startswith('NO_FILE'):
        p = Path(args.rdrpcatch)
        if p.exists() and p.stat().st_size > 0:
            print("Loading RdRPCATCH...", file=sys.stderr)
            rdrpcatch_df = load_rdrpcatch(args.rdrpcatch)

    rnavirhost_df = None
    if args.rnavirhost and not args.rnavirhost.startswith('NO_FILE'):
        p = Path(args.rnavirhost)
        if p.exists() and p.stat().st_size > 0:
            print("Loading RNAVirHost...", file=sys.stderr)
            rnavirhost_df = load_rnavirhost(args.rnavirhost)

    # --- Load BLASTn results ---
    blastn_dbs = {}
    if args.blastn_dir and not args.blastn_dir.startswith('NO_FILE'):
        p = Path(args.blastn_dir)
        if p.exists():
            print("Loading BLASTn results...", file=sys.stderr)
            blastn_dbs = load_all_blastn(args.blastn_dir)
            print(f"  Found {len(blastn_dbs)} database(s): "
                  f"{sorted(blastn_dbs.keys())}", file=sys.stderr)

    # --- Entrez lookup for core_nt ---
    entrez_df = None
    if 'core_nt' in blastn_dbs:
        print("Running Entrez host/lineage lookup for core_nt hits...",
              file=sys.stderr)
        # Determine the tname column for core_nt (e.g. tname_core_nt)
        corent_df = blastn_dbs['core_nt']
        tname_col = next(
            (c for c in corent_df.columns if c.startswith('tname_')), 'tname_nt'
        )
        entrez_df = get_host_lineage(
            corent_df,
            cache_path=args.host_lineage_cache,
            target_col=tname_col
        )

    # --- Report and merge ---
    report(genomad_df, checkv_df, rdrpcatch_df, rnavirhost_df)

    print("Merging annotations...", file=sys.stderr)
    merged_df = merge_master_table(
        checkv_df=checkv_df,
        genomad_df=genomad_df,
        rdrpcatch_df=rdrpcatch_df,
        blastn_dbs=blastn_dbs,
        entrez_df=entrez_df,
        rnavirhost_df=rnavirhost_df,
        ictv_fam_df=ictv_fam_df,
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(out, sep='\t', index=False)
    print(f"[merge_annotations] {len(merged_df)} rows → {out}", file=sys.stderr)


if __name__ == '__main__':
    main()
