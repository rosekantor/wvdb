#!/bin/bash
#SBATCH -t 00-72:00:00
#SBATCH -N 1
#SBATCH -o stdout_genomad_%j
#SBATCH -e stderr_genomad_%j
#SBATCH --job-name genomad
#SBATCH --qos=exempt

source /usr/workspace/kantor4/miniconda3/etc/profile.d/conda.sh
conda activate /usr/WS1/mlbiomon/software/mamba/envs/genomad

genomad_db=/p/vast1/mlbiomon/ref_data/genomad_db_1.9/genomad_db
in=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/votus_final.fasta
out=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/genomad/
genomad end-to-end --cleanup $in $out $genomad_db -t 120
