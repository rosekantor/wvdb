#!/usr/bin/env python3
'''
build taxonomy table for WVDB
'''
import argparse
import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
from pathlib import Path


def check_files(workdir, ictv_fam_file):
    """
    Validate required input files and raise one combined error if any are missing,
    not files, or empty.

    Returns
    -------
    tuple[pathlib.Path, ...]
        Paths to all validated files.
    """
    workdir = Path(workdir)

    errors = []

    if not workdir.exists():
        errors.append(f"Workdir does not exist: {workdir}")
    elif not workdir.is_dir():
        errors.append(f"Workdir is not a directory: {workdir}")

    checkv = workdir / "votus_final_quality_summary.tsv"
    genomad = workdir / "votus_final_virus_summary_genomad.tsv"
    rdrpcatch = workdir / "votus_final_rdrpcatch_output_annotated.tsv"
    core_nt_out = workdir / "votus_final-vs-corent.ani.tsv"
    imgvr_out = workdir / "votus_final-vs-IMGVR.ani.tsv"
    rnavirhost_csv = workdir / "votus_final_rnavirhost.csv"
    ictv_fam = Path(ictv_fam_file)

    required_files = {
        "checkv": checkv,
        "genomad": genomad,
        "rdrpcatch": rdrpcatch,
        "core_nt_out": core_nt_out,
        "imgvr_out": imgvr_out,
        "rnavirhost_csv": rnavirhost_csv,
        "ictv_fam": ictv_fam,
    }

    for name, path in required_files.items():
        if not path.exists():
            errors.append(f"{name}: file not found: {path}")
        elif not path.is_file():
            errors.append(f"{name}: expected a file, found something else: {path}")
        elif path.stat().st_size <= 0:
            errors.append(f"{name}: file is empty: {path}")

    if errors:
        raise FileNotFoundError("Input validation failed:\n  - " + "\n  - ".join(errors))

    return checkv, genomad, rdrpcatch, core_nt_out, imgvr_out, rnavirhost_csv, ictv_fam


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
    return genomad_tax_df


def load_rdrpcatch(rdrpcatch_file):
    rdrp_df = pd.read_csv(rdrpcatch_file, sep='\t')
    rdrp_df = rdrp_df.rename(columns={'Contig_name': 'contig'})
    rdrp_df['MMseqs_Taxonomy_2bLCA'] = (
        rdrp_df['MMseqs_Taxonomy_2bLCA'].fillna('').str.replace(' ', '_', regex=False)
    )

    tax = rdrp_df['MMseqs_Taxonomy_2bLCA']
    rdrp_df['Domain'] = tax.str.extract(r'd_([^;]+);?')
    rdrp_df['Realm'] = tax.str.extract(r'd_[^;]+;-_([^;]+);?')
    rdrp_df['Kingdom'] = tax.str.extract(r';k_([^;]+);?')
    rdrp_df['Phylum'] = tax.str.extract(r';p_([^;]+);?')
    rdrp_df['Class'] = tax.str.extract(r';c_([^;]+);?')
    rdrp_df['Order'] = tax.str.extract(r';o_([^;]+);?')
    rdrp_df['Family'] = tax.str.extract(r';f_([^;]+);?')
    rdrp_df['Genus'] = tax.str.extract(r';g_([^;]+);?')
    rdrp_df['Species'] = tax.str.extract(r';s_([^;]+);?')
    rdrp_df['Lowest'] = tax.str.extract(r'_([^;]+)$')

    rdrp_tax_df = rdrp_df[
        ['contig', 'Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order',
         'Family', 'Genus', 'Species', 'Lowest', 'MMseqs_Taxonomy_2bLCA', 'Best_hit_bitscore']
    ].copy()

    return rdrp_tax_df


def load_rnavirhost(rnavirhost_csv):
    rnavirhost_df = pd.read_csv(rnavirhost_csv, skiprows=1, names=['contig', 'Order', 'host_domain', 'host_chordata', 'evidence'])
    rnavirhost_df = rnavirhost_df.replace(np.nan, 'unknown')
    rnavirhost_df.rename(columns={'host_domain':'rnavirhost_domain', 
                                  'host_chordata':'rnavirhost_chordata'}, inplace=True)
    rnavirhost_df = rnavirhost_df[['contig', 
                                   'rnavirhost_domain', 
                                   'rnavirhost_chordata']]
    return rnavirhost_df
    

def process_blastn(tsv: str, db: str):
    '''
    starting with blast result ani file output from nayfach script (tsv) against a specific db (name of db)
    determine whether hits were genus or species level
    if there were multiple hits, select just the best one
    '''
    df = pd.read_csv(tsv, sep='\t')

    # filter for high identity and high coverage
    df['hit_level'] = 'none' 
    df.loc[(df.pid >=70) & (df.qcov >=85), 'hit_level'] = 'genus'
    df.loc[(df.pid >=95) & (df.qcov >=85), 'hit_level'] = 'species'
    df['idXcov'] = df.pid * df.qcov
    # if there are multiple hits, take the best one (groupby and idxmax)
    idx_max_values = df.groupby('qname')['idXcov'].idxmax()
    df = df.loc[idx_max_values]
    df = df.rename(columns={'qname':'contig',
                            'tname':f'tname_{db}',
                            'pid':f'pid_{db}',
                            'qcov':f'qcov_{db}',
                            'hit_level': f'hit_level_{db}'})
    df = df.drop(columns=['num_alns','tcov', 'idXcov'])

    return df


def batch_extract_hosts_lineages(genbank_ids):
    '''for a list of genbank ids, use a batch Entrez job to get 
    taxonomy, host, and isolation source if available'''
    
    ids_str = ",".join(genbank_ids)
    handle = Entrez.efetch(db="nucleotide", id=ids_str, rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    results = []
    for record in records:
        # Extract hosts and isolation source
        hosts = []
        isolation_source = None
        for feature in record.features:
            # note: not all records have host and isolation source
            if feature.type == "source":
                if "host" in feature.qualifiers:
                    hosts.extend(feature.qualifiers["host"])
                if "isolation_source" in feature.qualifiers:
                    isolation_source = feature.qualifiers.get("isolation_source", [None])[0]
        
        hosts = ','.join(hosts) if hosts else None
        isolation_source = isolation_source if isolation_source else None
        
        # Extract lineage
        lineage = record.annotations.get("taxonomy", None)
        lineage = ';'.join(lineage) # make a string from the list
        # combine all and save
        results.append((record.id, hosts, isolation_source, lineage))
    return results


def get_host_lineage_blasthit(corent, target_col='tname', batch_size=100):
    '''
    get host and lineage of best blast hit, starting with an ani tsv from blastn output
    ani tsv generated by nayfach ani script
    return an updated ani tsv with host, isolation source, and lineage of the hit as listed in the genbank entry
    also write these new columns to their own tsv
    '''
    genbank_id_list = list(corent[target_col].unique()) # get dereplicated list of genbank_ids hit by blast nt
    blastn_results = []
    for i in range(0, len(genbank_id_list), batch_size): # NCBI recommends <=500 per request
        print(f'querying entrez with batch {i}')
        
        batch_ids = genbank_id_list[i:i+batch_size]
        blastn_results.extend(batch_extract_hosts_lineages(batch_ids))
    
    # Convert results to DataFrame and merge with original
    blastn_results_df = pd.DataFrame(blastn_results, columns=[target_col, 'hosts_ntBlastHit', 'isolation_source_ntBlastHit', 'lineage_ntBlastHit'])
    blastn_results_df['lowest_tax_blasthit'] = blastn_results_df.lineage_ntBlastHit.str.extract(r'.*;(.*$)')
    blastn_results_df = corent[['contig', target_col]].merge(blastn_results_df, on=target_col, how='left')
    
    return blastn_results_df


def load_tax_tables(checkv, genomad, rdrpcatch, core_nt_out, imgvr_out, host_lineage_out, rnavirhost_csv, ictv_fam):
    
    host_lineage_out = Path(host_lineage_out)

    checkv_df = load_checkv(checkv)
    genomad_tax_df = load_genomad(genomad)
    rdrp_tax_df = load_rdrpcatch(rdrpcatch)
    rnavirhost_df = load_rnavirhost(rnavirhost_csv)

    corent = process_blastn(core_nt_out, 'nt')
    imgvr = process_blastn(imgvr_out, 'imgvr')

    for_entrez = corent[corent['hit_level_nt'].isin(['genus', 'species'])].copy()

    if host_lineage_out.is_file() and host_lineage_out.stat().st_size > 0:
        blastn_results_df = pd.read_csv(host_lineage_out, sep='\t')
    else:
        if len(for_entrez) == 0:
            blastn_results_df = pd.DataFrame(columns=['contig', 'tname_nt'])
        else:
            blastn_results_df = get_host_lineage_blasthit(for_entrez, target_col='tname_nt')

        host_lineage_out.parent.mkdir(parents=True, exist_ok=True)
        blastn_results_df.to_csv(host_lineage_out, sep='\t', index=False)

    ictv_fam_df = pd.read_csv(ictv_fam, sep='\t')

    return checkv_df, genomad_tax_df, rdrp_tax_df, corent, imgvr, blastn_results_df, rnavirhost_df, ictv_fam_df


def merge_master_table(
    checkv_df,
    genomad_tax_df,
    rdrp_tax_df,
    corent,
    imgvr,
    blastn_results_df,
    ictv_fam_df,
    rnavirhost_df,
):
    """
    Merge all annotation sources into a final master table.

    Returns
    -------
    pd.DataFrame
    """
    genomad_merge = genomad_tax_df[
        ['contig', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']
    ].drop_duplicates('contig')

    rdrp_merge = rdrp_tax_df[
        ['contig', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Lowest']
    ].drop_duplicates('contig')

    checkv_merge = checkv_df.drop_duplicates('contig')
    rnavirhost_merge = rnavirhost_df.drop_duplicates('contig')

    taxhost_df = pd.merge(
        genomad_merge,
        rdrp_merge,
        on='contig',
        how='outer',
        suffixes=('_gNd', '_RdRp')
    )

    taxhost_df = taxhost_df.rename(
        columns={'Genus': 'Genus_RdRp', 'Species': 'Species_RdRp', 'Lowest': 'Lowest_RdRp'}
    )

    taxhost_df = taxhost_df.merge(checkv_merge, on='contig', how='outer')
    taxhost_df = taxhost_df.merge(rnavirhost_merge, on='contig', how='outer')
    taxhost_df['rnavirhost_domain'] = taxhost_df['rnavirhost_domain'].fillna('Unknown')

    taxhost_df = taxhost_df.merge(corent, on='contig', how='left')
    taxhost_df = taxhost_df.merge(blastn_results_df, on=['contig', 'tname_nt'], how='left')
    taxhost_df = taxhost_df.merge(imgvr, on='contig', how='left')

    taxhost_df['Order_consensus'] = taxhost_df['Order_RdRp']
    taxhost_df.loc[taxhost_df['Order_consensus'].isna(), 'Order_consensus'] = taxhost_df['Order_gNd']
    taxhost_df['Order_consensus'] = taxhost_df['Order_consensus'].fillna('Unclassified')

    taxhost_df['Family_consensus'] = taxhost_df['Family_RdRp']
    taxhost_df.loc[taxhost_df['Family_consensus'].isna(), 'Family_consensus'] = taxhost_df['Family_gNd']
    taxhost_df['Family_consensus'] = taxhost_df['Family_consensus'].fillna('Unclassified')

    taxhost_df = taxhost_df.merge(
        ictv_fam_df,
        how='left',
        left_on='Family_consensus',
        right_on='Family_ICTV'
    )
    taxhost_df['host_ICTV_simple'] = taxhost_df['host_ICTV_simple'].fillna('Unknown')

    return taxhost_df


def make_master_table(
    checkv_df,
    genomad_tax_df,
    rdrp_tax_df,
    corent,
    imgvr,
    blastn_results_df,
    ictv_fam_df,
    rnavirhost_df,
    outfile,
):
    """
    Orchestrate reporting, RNAVirHost preparation, merging, and output writing.

    Returns
    -------
    pd.DataFrame
    """
    provirus_df = genomad_tax_df[genomad_tax_df['contig'].str.contains('provirus', na=False)]
    high_kmer_df = checkv_df[checkv_df['kmer_freq'] > 1.2]

    if len(provirus_df) > 0:
        print(
            f"geNomad reports {len(provirus_df)} proviruses. "
            f"Consider trimming before continuing: {provirus_df['contig'].tolist()}"
        )

    if len(high_kmer_df) > 0:
        print(
            f"checkV reports {len(high_kmer_df)} contigs with kmer_freq > 1.2. "
            f"Consider trimming before continuing: {high_kmer_df['contig'].tolist()}"
        )

    if len(rdrp_tax_df) != len(rdrp_tax_df['contig'].unique()):
        print("RdRpCATCH output has multiple RdRp genes per contig")

    print(
        f"Annotation counts\n"
        f"checkV: {len(checkv_df)}\n"
        f"RdRpCATCH: {len(rdrp_tax_df)}\n"
        f"geNomad: {len(genomad_tax_df)}"
    )

    print(f"RNAVirHost: {len(rnavirhost_df)}")

    taxhost_df = merge_master_table(
        checkv_df=checkv_df,
        genomad_tax_df=genomad_tax_df,
        rdrp_tax_df=rdrp_tax_df,
        corent=corent,
        imgvr=imgvr,
        blastn_results_df=blastn_results_df,
        ictv_fam_df=ictv_fam_df,
        rnavirhost_df=rnavirhost_df,
    )

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    taxhost_df.to_csv(outfile, sep='\t', index=False)

    return taxhost_df


def main():
    parser = argparse.ArgumentParser(
        description="Build a merged viral annotation table from checkV, geNomad, RdRpCATCH, ANI, ICTV, and RNAVirHost outputs."
    )
    parser.add_argument(
        "workdir",
        type=Path,
        help="Working directory containing expected input files, and where outputs will be written.",
    )
    parser.add_argument(
        "--ictv-fam",
        type=Path,
        required=True,
        help="Path to ICTV family metadata TSV.",
    )
    parser.add_argument(
        "--entrez-email",
        required=True,
        help="Email address to use for NCBI Entrez requests.",
    )

    args = parser.parse_args()

    workdir = args.workdir
    ictv_fam = args.ictv_fam
    entrez_email = args.entrez_email

    # Set Entrez email for downstream NCBI lookups
    Entrez.email = entrez_email

    checkv, genomad, rdrpcatch, core_nt_out, imgvr_out, rnavirhost_csv, ictv_fam = check_files(workdir, ictv_fam)

    host_lineage_out = workdir / "votus_final_host_lineage.tsv"
    outfile = workdir / "votus_final_full_annotation.tsv"

    (
        checkv_df,
        genomad_tax_df,
        rdrp_tax_df,
        corent,
        imgvr,
        blastn_results_df,
        rnavirhost_df,
        ictv_fam_df,
    ) = load_tax_tables(
        checkv=checkv,
        genomad=genomad,
        rdrpcatch=rdrpcatch,
        core_nt_out=core_nt_out,
        imgvr_out=imgvr_out,
        host_lineage_out=host_lineage_out,
        rnavirhost_csv=rnavirhost_csv,
        ictv_fam=ictv_fam,
    )

    make_master_table(
        checkv_df=checkv_df,
        genomad_tax_df=genomad_tax_df,
        rdrp_tax_df=rdrp_tax_df,
        corent=corent,
        imgvr=imgvr,
        blastn_results_df=blastn_results_df,
        ictv_fam_df=ictv_fam_df,
        rnavirhost_df=rnavirhost_df,
        outfile=outfile,
    )


if __name__ == "__main__":
    main()

# workdir = "/Users/kantor4/data/trilab_ldrd/virus_db_analysis/clusters_v3"

# # Set email (required by NCBI)
# Entrez.email = "kantor4@llnl.gov"


