# wvdb

## Description
This project contains two Nextflow pipelines for building and annotating virus genome databases from metagenomics: `wvdb-build` and `wvdb-annotate`. These workflows were developed for the Wastewater Virus Database project (WVDB).

`wvdb-build` is designed to take as input a set of prefiltered, high-quality viral genomes recovered from across multiple metagenomic samples (locations, timepoints). Clustering and trimming is used in an attempt to remove potential chimeras resulting from misassemblies. It produces a dereplicated, trimmed set of unique vOTU-level genomes to be used for downstream analyses.     

`wvdb-annotate` takes the final vOTU fasta file and annotates the viral genomes with respect to taxonomy and predicted host, as well as matches to a user-defined set of databases. It produces a flat annotation table to be used in downstream analyses.    

For more information on installation and usage, please see READMEs within each workflow.

## Authors and acknowledgment
This code was developed by Rose Kantor at Lawrence Livermore National Laboratory, with input from Migun Shakya (LANL), and Joe Wakim (LLNL). Code and READMEs were written with assistance from Claude Sonnet 5.    
This code also makes use of a BLAST-based method for average nucleotide identity, developed by Nayfach et al. as part of UHGV (https://github.com/snayfach/MGV/blob/master/ani_cluster/blastani.py).
