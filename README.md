# wvdb

## Description
This project contains two Nextflow pipelines for building and annotating virus genome databases from metagenomics: `wvdb-build` and `wvdb-annotate`. These workflows were developed for the Wastewater Virus Database project (WVDB).

`wvdb-build` is designed to take as input a set of prefiltered, high-quality viral genomes recovered from across multiple metagenomic samples (locations, timepoints). Clustering and trimming is used in an attempt to remove potential chimeras resulting from misassemblies. It produces a dereplicated, trimmed set of unique vOTU-level genomes to be used for downstream analyses.     

`wvdb-annotate` takes the final vOTU fasta file and annotates the viral genomes with respect to taxonomy and predicted host, as well as matches to a user-defined set of databases. It produces a flat annotation table to be used in downstream analyses.    

For more information on installation and usage, please see READMEs within each workflow.

## Authors and acknowledgment
This code was developed by Rose Kantor at Lawrence Livermore National Laboratory, with input from Migun Shakya (LANL), and Joe Wakim (LLNL). Code and READMEs were written with assistance from Claude Sonnet 5.    
This code also makes use of a BLAST-based method for average nucleotide identity, developed by Nayfach et al. as part of UHGV (https://github.com/snayfach/MGV/blob/master/ani_cluster/blastani.py).

## Release
This code was released under LLNL release number LLNL-CODE-2024371.

## Notice
This work was produced under the auspices of the U.S. Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.

Neither the United States Government nor Lawrence Livermore National Security, LLC, nor any of their employees makes any warranty, expressed or implied, or assumes any legal liability or responsibility for the
accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately owned
rights. Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, LLC. The views and
opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence Livermore National Security, LLC,
and shall not be used for advertising or product endorsement purposes.
