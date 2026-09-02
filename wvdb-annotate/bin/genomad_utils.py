#!/usr/bin/env python3
"""
genomad_utils.py — shared helpers for parsing geNomad output.

geNomad appends '|provirus_<start>_<end>' to seq_name in its
<prefix>_virus_summary.tsv for any sequence where it identified an
integrated provirus region within a larger contig — one row per contig,
never both a plain-ID row and a provirus-suffixed row for the same
contig. This suffixed ID does not match the original FASTA header,
which breaks any downstream join keyed on exact contig ID. In practice
this surfaced as duplicate rows in RNAVirHost's orders.csv: RdRPCATCH
and the FASTA itself use the plain contig ID, while geNomad's row for
the same contig used the suffixed ID, so an outer join produced two
separate rows (one matching RdRPCATCH/FASTA, one matching only geNomad)
where there should have been exactly one.

CheckV does not always catch these same regions — a genome can carry
clear proviral/host gene content (per geNomad) while still passing
CheckV's own contamination checks, especially when the same integrated
element is reproducibly assembled across multiple independent samples
(suggesting a real phage-plasmid or temperate phage with genuine host
genes, not a chimeric assembly artifact). These cases are flagged, not
discarded — the geNomad classification is preserved and the original
FASTA sequence is kept as-is.

This module normalizes the suffix away for joining purposes, while
preserving the provirus flag and coordinates for annotation.
"""

import re

import numpy as np
import pandas as pd


_PROVIRUS_SUFFIX_RE = re.compile(r'^(.+)\|provirus_(\d+)_(\d+)$')


def normalize_genomad_id(seq_name):
    """
    Split a geNomad seq_name into (base_contig_id, is_provirus, start, end).

    '<id>|provirus_<start>_<end>' → ('<id>', True, start, end)
    '<id>'                        → ('<id>', False, None, None)
    """
    m = _PROVIRUS_SUFFIX_RE.match(str(seq_name))
    if m:
        return m.group(1), True, int(m.group(2)), int(m.group(3))
    return seq_name, False, None, None


def load_genomad_virus_summary(path):
    """
    Load geNomad's <prefix>_virus_summary.tsv and normalize contig IDs so
    they match the original FASTA headers exactly, regardless of whether
    geNomad reported a plain contig or a '|provirus_<start>_<end>' region.
    geNomad emits exactly one row per contig, so no deduplication is
    needed here — normalization alone resolves the ID mismatch.

    Adds columns:
      contig                  — normalized ID, matches FASTA header
      genomad_provirus        — bool, True if this row is a provirus region
      genomad_provirus_start  — int or NaN
      genomad_provirus_end    — int or NaN
    """
    df = pd.read_csv(path, sep='\t')

    normalized = df['seq_name'].apply(normalize_genomad_id)
    df['contig']                 = [n[0] for n in normalized]
    df['genomad_provirus']       = [n[1] for n in normalized]
    df['genomad_provirus_start'] = [n[2] for n in normalized]
    df['genomad_provirus_end']   = [n[3] for n in normalized]

    # Split taxonomy string into rank columns
    ranks = ['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']
    df[ranks] = df['taxonomy'].str.split(';', expand=True)
    df[ranks] = df[ranks].replace('', np.nan)

    return df
