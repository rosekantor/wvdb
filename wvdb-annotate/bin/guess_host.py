import pandas as pd
import numpy as np

# define inputs and outputs
taxhost_file = '/Users/kantor4/data/trilab_ldrd/virus_db_analysis/clusters_v3/votus_final_full_annotation.tsv' # from merge_annotations.py
final_table = '/Users/kantor4/data/trilab_ldrd/virus_db_analysis/clusters_v3/votus_final_full_annotation_hostupdate.tsv'

taxhost_df = pd.read_csv(taxhost_file, sep='\t')

taxhost_df.loc[taxhost_df.host_ICTV.isna(), 'host_ICTV'] = 'unknown'

# make lowercase
taxhost_df.rnavirhost_domain = taxhost_df.rnavirhost_domain.str.lower()
taxhost_df.host_ICTV = taxhost_df.host_ICTV.str.lower()

# adjust host names to match between tools
taxhost_df.loc[taxhost_df.rnavirhost_domain == 'viridiplantae', 'rnavirhost_domain'] = 'plants'
taxhost_df.loc[taxhost_df.rnavirhost_domain == 'chordata', 'rnavirhost_domain'] = 'vertebrates'
taxhost_df.loc[taxhost_df.rnavirhost_domain == 'invertebrate', 'rnavirhost_domain'] = 'invertebrates'

# Get info from core-nt blast results where possible
# Normalize BLASTn core-nt host results (pulled down with Entrez during merge_annotations.py)
# 1. use LLM to check if isolation source is fecal related
#   a. feed the following list to LLM: list(taxhost_df[taxhost_df.hosts_ntBlastHit != ''].isolation_source_ntBlastHit.unique())
#   b. request isolation_category=['fecal associated', 'digestive tract non-fecal', 'plant-associated', 'other']

# example of read-in dict
# isolation_source_dict = {
#     'human feces sample in Mie, Japan': 'fecal associated',
#     'fecal samples': 'fecal associated',
#     'stool': 'fecal associated',
#     np.nan: 'other',
#     'freshwater mussel tissue biopsies': 'other',
#     'rice field': 'plant-associated',
#     'Rice Field': 'plant-associated',
#     'sewage': 'other',
#     'fecal sample': 'fecal associated',

# 2. use LLM to categorize host, using API key
#   a. feed the following list to LLM: list(taxhost_df[taxhost_df.hosts_ntBlastHit != ''].hosts_ntBlastHit.unique())
#   b. request host_categories = ['bacteria', 'unknown', 'vertebrates', 'plants', 'invertebrates', 'fungi']
#   c. bring back a simplified dict for manual review (at least give option for this)

# example of read-in dict
# blastnt_host_dict = {
#     'Homo sapiens': 'vertebrates',
#     'Sus scrofa': 'vertebrates',
#     'pepper': 'plants',
#     'Enterobacteriaceae': 'bacteria',
#     'bird': 'vertebrates',
#     'Pseudomonas putida KT2440 + pKJK5': 'bacteria',
#     'Channeled applesnail': 'invertebrates',
#     'Ortmanniana pectorosa': 'invertebrates',...}

isolation_source_df = pd.DataFrame(
    list(isolation_source_dict.items()),
    columns=["isolation_source_ntBlastHit", "isolation_source_ntBlastHit_simple"]
)
blastnt_host_df = pd.DataFrame(
    list(blastnt_host_dict.items()),
    columns=["hosts_ntBlastHit", "hosts_ntBlastHit_simple"]
)
# apply adjustments as manually determined
blastnt_host_df[blastnt_host_df.hosts_ntBlastHit_simple=='unknown']
# 'Guinardia delicatula' is a diatom --> 'invertebrate'
# 'crida cinerea' should be Acrida --> 'invertebrate'
blastnt_host_df.loc[blastnt_host_df.hosts_ntBlastHit == 'Guinardia delicatula', 'hosts_ntBlastHit_simple'] = 'invertebrates'
blastnt_host_df.loc[blastnt_host_df.hosts_ntBlastHit == 'crida cinerea', 'hosts_ntBlastHit_simple'] = 'invertebrates'

# add back to tax_host_df
taxhost_df = taxhost_df.merge(isolation_source_df, how='left', on='isolation_source_ntBlastHit')
taxhost_df = taxhost_df.merge(blastnt_host_df, how='left', on='hosts_ntBlastHit')

# apply decision tree logic
# * compare ICTV vs RNAVirHost, determine agreement
# * where BLAST exists:
#     * for ICTV-RNAVirHost agree, check against BLAST --> if disagree, review
#     * for ICTV-RNAVirHost disagree, check against BLAST --> if BLAST agrees with one, use BLAST
#     * for ICTV-RNAVirHost unknown, choose BLAST but flag if isolation source is fecal (could be dietary)
# * where no BLAST exists:
#     * for ICTV-RNAVirHost agree --> use RNAVirHost
#     * for ICTV-RNAVirHost disagree and ICTV known --> use ICTV
#     * for ICTV-RNAVirHost disagree and RNAVirHost known --> use RNAVirHost


# Assumptions:
# - df["rnavirhost_domain"] uses "unknown" for unknown values
# - df["host_ICTV"] uses "unknown" for unknown values
# - df["hosts_ntBlastHit_simple"] uses np.nan for unknown values
# - df["isolation_source_ntBlastHit_simple"] may contain strings such as "fecal associated"

def ictv_vs_rnavirhost_agreement(row):
    """
    Compare RNAVirHost to ICTV.
    ICTV may contain multiple comma-separated values.
    """
    rna = str(row["rnavirhost_domain"]).strip()
    ictv = str(row["host_ICTV"]).strip()

    if rna == "unknown" or ictv == "unknown":
        return "unknown"

    ictv_vals = [x.strip() for x in ictv.split(",")]
    return "agree" if rna in ictv_vals else "disagree"


def blast_exists(row):
    return pd.notna(row["hosts_ntBlastHit_simple"])


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

    # first split: do RNAVirHost and ICTV method have known values and agree?
    if ictv_rna_agree:
        # 1) if yes, final host = RNAVirHost host, method = ICTV and RNAVirHost
        # 1a) if BLAST hit exists and agrees, add BLAST to the method list
        if blast_agrees_ictv:
            return out(blast_host, "ICTV_RNAVirHost_ntBlastHit_host")

        # 1a) if BLAST hit exists and disagrees, flag for review
        if has_blast and not blast_agrees_ictv:
            return out(rna, "ICTV_RNAVirHost", True, "ICTV and RNAVirHost agree, BLAST disagrees")

        return out(rna, "ICTV_RNAVirHost")

    # 2) if both are known but disagree, check the following:
    if both_known and not ictv_rna_agree:
        # prioritize agreement between ICTV and BLAST if present
        if has_blast and blast_agrees_ictv:
            return out(blast_host, "ICTV_ntBlastHit_host")

        # 2a) if BLAST hit doesn't exist, final host = ICTV method
        if not has_blast:
            return out(ictv, "ICTV")

        # 2b) if BLAST hit exists, check whether BLAST hit agrees with ICTV or RNAVirhost
        # 2b1) if BLAST hit does agree with one of these, use the BLAST hit, note both agreeing methods
        if blast_agrees_rna:
            return out(blast_host, "RNAVirHost_ntBlastHit_host")

        # 2b2) if BLAST hit disagrees with both ICTV and RNAVirHost, use ICTV method and flag for review
        return out(ictv, "ICTV", True, "ICTV and RNAVirHost disagree, BLAST disagrees with both")

    # 3) if either RNAVirHost or ICTV method is unknown:
    if not both_known:
        # 3c) if both are unknown and BLAST hit exists, use BLAST hit
        # isolation source is checked only when BLAST is used alone
        if both_unknown and has_blast:
            if str(blast_iso) == "fecal associated":
                return out(blast_host, "ntBlastHit_host", True, "BLAST used alone, isolation source is fecal associated")
            return out(blast_host, "ntBlastHit_host")

        # 3d) if both are unknown and BLAST hit does not exist, final host is unknown
        if both_unknown and not has_blast:
            return out("unknown", "none", False, None)

        # 3a) if BLAST hit exists, check for agreement with the known (RNAVirHost or ICTV)
        if has_blast:
            # 3a1) if agreement, use BLAST hit, list all agreeing methods
            if ictv_known and blast_agrees_ictv:
                return out(blast_host, "ICTV_ntBlastHit_host")

            if rna_known and blast_agrees_rna:
                return out(blast_host, "RNAVirHost_ntBlastHit_host")

            # 3a2) if disagree and ICTV is known, use ICTV method
            if ictv_known:
                return out(ictv, "ICTV")

            # 3a3) if disagree and ICTV is unknown, use BLAST
            # isolation source is checked only when BLAST is used alone
            if not ictv_known:
                if str(blast_iso) == "fecal associated":
                    return out(blast_host, "ntBlastHit_host", True, "BLAST used alone, isolation source is fecal associated")
                return out(blast_host, "ntBlastHit_host")

        # 3b) if BLAST hit doesn't exist, use whichever is known (ICTV or RNAVirHost)
        if not has_blast:
            if ictv_known:
                return out(ictv, "ICTV")
            if rna_known:
                return out(rna, "RNAVirHost")

    return out("review", "review", True, "Unhandled case")


# Apply to dataframe
taxhost_df["ictv_rnavirhost_agreement"] = taxhost_df.apply(ictv_vs_rnavirhost_agreement, axis=1)
taxhost_df["blast_exists"] = taxhost_df.apply(blast_exists, axis=1)

taxhost_df[[
    "final_host_determination",
    "final_tool_used",
    "review_flag",
    "review_reason"
]] = taxhost_df.apply(choose_final, axis=1)

# add bacteria host for Leviviricetes where Order was unclassified
bact_classes = ['Leviviricetes', 'Caudoviricetes', 'Vidaverviricetes', 'Faserviricetes']
taxhost_df.loc[(taxhost_df.Order_consensus == 'Unclassified') & 
    ((taxhost_df.Class_gNd.isin(bact_classes)) | 
     (taxhost_df.Class_RdRp.isin(bact_classes))), 'final_host_determination'] = 'bacteria'

taxhost_df.loc[(taxhost_df.Order_consensus == 'Unclassified') & 
    ((taxhost_df.Class_gNd.isin(bact_classes)) | 
     (taxhost_df.Class_RdRp.isin(bact_classes))), 'final_tool_used'] = 'manual_by_class'

# review and fix obvious cases
# (Martelli are usually but not always plants, Norzi and Timlo are always bacteria)
taxhost_df[(taxhost_df.Order_consensus == 'Martellivirales') & (taxhost_df.final_host_determination != 'plants')]
## confirm if some are are actually not plants

# fix here
taxhost_df.loc[(taxhost_df.Order_consensus == 'Norzivirales') & (taxhost_df.final_host_determination != 'bacteria'), 'final_host_determination'] = 'bacteria'
taxhost_df.loc[(taxhost_df.Order_consensus == 'Timlovirales') & (taxhost_df.final_host_determination != 'bacteria'), 'final_host_determination'] = 'bacteria'

# Consolidate "other" as a host category
# Make an "other" category for hosts with very few vOTUs (less than 100)
# apply this category to the summary table and to the main taxhost_df table
summary_host = taxhost_df.groupby('final_host_determination')[['contig']].count().reset_index()
other = summary_host[summary_host.contig < 100].final_host_determination.to_list()

taxhost_df['host_final'] = taxhost_df['final_host_determination']
taxhost_df.loc[taxhost_df.final_host_determination.isin(other), 'host_final'] = 'other'

# write final table
taxhost_df.to_csv(final_table, sep='\t', index=False)