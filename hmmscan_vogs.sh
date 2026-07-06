#!/bin/bash
#SBATCH --job-name=hmmscan_vogs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=120
#SBATCH --time=2-00:00:00
#SBATCH --qos=exempt
#SBATCH --array=1-5
#SBATCH --output=stdout_hmmscan_%A_%a.log
#SBATCH --error=stderr_hmmscan_%A_%a.log

set -euo pipefail

indir="/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/prodigal/votus_final.genes.faa.split"
outdir="/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/hmm_out"
db="/p/vast1/mlbiomon/analysis/securebio/vogdb/vog234/vogs234.hmm"

mkdir -p "$outdir"

part=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")
input="${indir}/votus_final.genes.part_${part}.faa"
base="$(basename "$input" .faa)"

[[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 1; }

# Prefer pressing once before submitting the array.
# This check avoids rerunning if already present.
if [[ ! -f "${db}.h3m" || ! -f "${db}.h3i" || ! -f "${db}.h3f" || ! -f "${db}.h3p" ]]; then
    echo "Pressed HMM database files not found for: $db" >&2
    echo "Run: hmmpress \"$db\" before submitting this array." >&2
    exit 1
fi

hmmscan \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    --domtblout "$outdir/${base}-vs-vogs.domtblout.txt" \
    "$db" \
    "$input" \
    > /dev/null
