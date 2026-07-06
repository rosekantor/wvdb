#!/bin/bash
#SBATCH -t 00-24:00:00
#SBATCH -N 1
#SBATCH -o stdout_checkv_%j
#SBATCH -e stderr_checkv_%j
#SBATCH --job-name checkv

workdir=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3
checkv=/usr/WS1/mlbiomon/software/mamba/envs/checkv/bin/checkv
checkvdb=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/checkv_db/checkv-db-v1.5
$checkv end_to_end -d $checkvdb -t 100 --remove_tmp $workdir/votus_final.fasta $workdir/checkv
