#!/bin/bash
#SBATCH -t 00-72:00:00
#SBATCH --qos=exempt
#SBATCH -N 1
#SBATCH -o stdout_rdrpcatch_%j
#SBATCH -e stderr_rdrpcatch_%j
#SBATCH --job-name rdrpcatch

rdrpcatch=/usr/WS1/mlbiomon/software/rdrpcatch/bin/rdrpcatch
db=/p/vast1/mlbiomon/ref_data/rdrp_catch_db
fasta=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/votus_final.fasta
outdir=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/rdrpcatch/rdrpcatch_out

$rdrpcatch scan -i $fasta -o $outdir -db-dir $db --cpus 128 -seq_type nuc
