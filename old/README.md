# wf3-pipeline

Viral genome clustering and reference-guided trimming pipeline, converted from `run_wf3_v2.sh`.

Genomes from multiple samples are collected, clustered with vclust, and trimmed via three
complementary branches (cluster-pair trimming, second-rank trimming, and BLAST-mode trimming
against reference databases). The final output is a non-redundant set of viral genome
representatives produced by reclustering.

## Repository structure

```
wf3-pipeline/
├── main.nf                       # Entry point; wires all steps and subworkflows
├── nextflow.config               # Params, profiles (local, slurm, test)
├── conf/
│   └── base.config               # Per-process CPU/memory labels
├── modules/
│   ├── collect.nf                # Step 1: collect and merge input FASTAs
│   ├── vclust.nf                 # Steps 2 & 6: vclust prefilter/align/cluster + centroids
│   ├── trim_filter.nf            # Shared: parse_clusters, trim_genomes, checkv, completeness_filter
│   └── blast_trim.nf             # Branch C: blastn, blastani, trim_genomes (blast mode), merge
├── subworkflows/
│   ├── trim_and_filter.nf        # Branches A & B: shared trim→checkv→filter subworkflow
│   └── branchC.nf                # Branch C: parallel BLAST fan-out subworkflow
├── bin/                          # Python helper scripts (added to PATH by Nextflow automatically)
│   ├── blastani_nayfach.py
│   ├── parse_clusters_v2.py
│   ├── trim_genomes.py
│   └── completeness_filter.py
├── envs/
│   └── wf3.yml                   # Conda environment for all tools
└── test/
    └── data/                     # Small test inputs for local validation
```

## Dependencies

All tools are managed via the conda environment in `envs/wf3.yml`.

| Tool | Version | Notes |
|---|---|---|
| Nextflow | ≥ 23.10 | Requires Java 11+ |
| seqkit | 2.9.0 | |
| vclust | 1.3.1 | |
| blastn | 2.16.0+ | via `blast` conda package |
| checkv | latest (1.1.1+) | |
| nucmer | 4.0.1 | via `mummer4` conda package |
| diamond | ≥ 2.0.9 | checkv hard dependency |
| Python | 3.11 | biopython, pandas, numpy required by bin/ scripts |

## Installation

### 1. Create the conda environment

```bash
# Recommended: use mamba for faster dependency solving
conda install -n base -c conda-forge mamba

mamba env create -f envs/wf3.yml
conda activate wf3
```

### 2. Verify all tools resolve

Each tool uses a different flag for version output — use the commands below
rather than a generic loop:

```bash
seqkit version                    # seqkit v2.9.0
vclust -v                         # vclust 1.3.1
blastn -version                   # blastn: 2.16.0+
checkv -h 2>&1 | head -1          # checkv 1.1.x (no dedicated --version flag)
nucmer --version                  # 4.0.1
nextflow -v                       # nextflow version 23.x
```

Verify Python helper scripts are importable and executable from `bin/`:

```bash
for script in blastani_nayfach.py parse_clusters_v2.py trim_genomes.py completeness_filter.py; do
    python bin/$script --help > /dev/null 2>&1         && echo "OK: $script"         || echo "CHECK: $script — may not support --help; try running with no args"
done
```

### 3. Validate the pipeline DAG (no data needed)

```bash
nextflow run main.nf -stub -profile local
```

## Running the pipeline

### Local development / macOS

```bash
conda activate wf3

nextflow run main.nf -profile local \
    --fastqdir /path/to/fastq \
    --outdir   /path/to/results \
    --checkvdb /path/to/checkv-db-v1.5 \
    --refseq_ev_blastdb /path/to/esviritu_plus_refseq_virus.fna \
    --imgvr_blastdb     /path/to/IMGVR5_UViG.fna \
    --max_memory '32 GB'
```

### SLURM cluster (production)

```bash
conda activate wf3

nextflow run main.nf -profile slurm \
    --fastqdir /path/to/fastq \
    --outdir   /path/to/results \
    --checkvdb /path/to/checkv-db-v1.5 \
    --refseq_ev_blastdb /path/to/esviritu_plus_refseq_virus.fna \
    --imgvr_blastdb     /path/to/IMGVR5_UViG.fna \
    -resume
```

The `-resume` flag restarts from the last successful checkpoint if any process fails.

### Test profile (small data)

Place test FASTAs under `test/data/` following the expected layout:

```
test/data/
└── <sample>/
    └── assembly/
        └── <run>/
            ├── filtered_virus.fasta
            └── filtered_provirus.fasta
```

Then run:

```bash
nextflow run main.nf -profile test \
    --checkvdb /path/to/checkv-db-v1.5 \
    --refseq_ev_blastdb /path/to/esviritu_plus_refseq_virus.fna \
    --imgvr_blastdb     /path/to/IMGVR5_UViG.fna
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--fastqdir` | required | Top-level wf2 output dir containing `<sample>/assembly/<run>/filtered_*.fasta` |
| `--outdir` | required | Output directory for all results |
| `--checkvdb` | required | CheckV database directory (e.g. `checkv-db-v1.5`) |
| `--refseq_ev_blastdb` | required | RefSeq virus + EsViritu combined BLAST db (`.fna` + index files) |
| `--imgvr_blastdb` | required | IMG-VR v5 (MetaVR) BLAST db (`.fna` + index files) |
| `--ani` | 0.95 | ANI threshold for vclust clustering |
| `--qcov` | 0.85 | Query coverage threshold for vclust |
| `--completeness` | 90 | CheckV completeness threshold (%) |
| `--threads` | 100 | Thread count for high-CPU processes (overridden by profiles) |
| `--max_memory` | `'64 GB'` | Memory cap for high-CPU processes; set lower for local runs |

## Output structure

```
outdir/
├── 1_filtered/
│   ├── filtered_all.fasta              # all input genomes merged
│   ├── filtered_all.length.txt         # per-sequence lengths
│   ├── filtered_virus_all.fasta
│   └── filtered_provirus_all.fasta
├── 2_clustered/
│   ├── vclust_clusters.tsv             # cluster assignments
│   └── vclust_ani.ids.tsv
├── 3_branchA/                          # rank-1 vs rank-2 cluster trimming
│   ├── cluster_pairs.tsv
│   ├── trimmed.fasta
│   ├── checkv_out/
│   ├── cluster_reps_complete.fasta     # → fed to step 6
│   └── cluster_reps_incomplete_or_na.txt
├── 4_branchB/                          # rank-2 vs rank-3 trimming (A incompletes only)
│   ├── cluster_pairs.tsv
│   ├── trimmed.fasta
│   ├── checkv_out/
│   ├── cluster_reps23_complete.fasta   # → fed to step 6
│   └── cluster_reps_incomplete.fasta   # → fed to branch C
├── 5_branchC/                          # BLAST-mode trimming (singletons + B incompletes)
│   ├── branchC_genomes.fasta
│   ├── refseq_ev/                      # blastn + blastani + trim outputs
│   ├── metavr/
│   ├── branchC_trimmed.fasta
│   ├── checkv_out/
│   └── blastmode_complete.fasta        # → fed to step 6
└── 6_reclustered/
    ├── recluster_input.fasta           # merged A + B + C complete sets
    ├── vclust_clusters.tsv
    ├── vclust_centroids.txt
    └── vclust_centroids.fasta          # ← FINAL OUTPUT
```

## Pipeline overview

```
Input FASTAs (virus + provirus per sample)
        │
        ▼
   Step 1: Collect & merge
        │
        ▼
   Step 2: vclust cluster (prefilter → align → cluster)
        │
   ┌────┴──────────────────────────┐
   ▼                               ▼
Branch A (rank 1→2)         Branch B (rank 2→3)        Branch C (singletons + B incompletes)
trim → checkv → filter      trim → checkv → filter     BLAST ×2 → trim ×2 → merge → checkv → filter
   │                               │                               │
   └───────────────────────────────┴───────────────────────────────┘
                                   │
                                   ▼
                          Step 6: Recluster (vclust)
                                   │
                                   ▼
                          vclust_centroids.fasta
```

## Development notes

### Conda environment updates

To add or update a package:
1. Edit `envs/wf3.yml`
2. Run `mamba env update -f envs/wf3.yml --prune`
3. Re-run stub validation: `nextflow run main.nf -stub -profile local`
4. Commit the updated yml

### Platform notes

- All tools support **linux-64** (cluster) and **macOS ARM64** (Apple Silicon).
- `vclust=1.3.1` supports macOS ARM64 natively — no Rosetta needed.
- `mummer4` is the bioconda package name for nucmer v4; the binary is still called `nucmer`.
- The `local` profile caps CPUs and memory to what the machine actually has, using
  `Runtime.runtime.availableProcessors()` and `--max_memory`.

### publishDir and task.ext

Intermediate outputs are published via `publishDir { "${params.outdir}/${task.ext.publish_dir}" }`.
The `task.ext.publish_dir` values are set in `nextflow.config` using `withName` selectors scoped
to each subworkflow alias (`BRANCH_A:*`, `BRANCH_B:*`, `BRANCH_C:*`). The closure form is required
because `task.ext` is only accessible at runtime, not at compile time.
