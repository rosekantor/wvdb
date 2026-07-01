# wvdb_build

Viral genome database build pipeline. Collects viral and proviral genome
assemblies across samples, clusters them with vclust, trims chimeric
assemblies using pairwise alignments, assesses completeness with CheckV,
and produces a non-redundant set of representative genomes via reclustering.

This pipeline lives in the `wvdb-build/` subdirectory of the wvdb repository,
alongside `wvdb-annotate/` and other workflows.

---

## Repository structure

```
wvdb-build/
├── main.nf                         # Entry point; wires all steps
├── nextflow.config                 # Params, profiles (local, slurm, conda, test)
├── conf/
│   └── base.config                 # Per-process CPU/memory resource labels
├── modules/
│   ├── collect.nf                  # Step 1: collect and merge input FASTAs
│   ├── vclust.nf                   # Steps 2 & 5: vclust prefilter/align/cluster + centroids
│   ├── cluster_trim_tools.nf       # Step 3 processes: parse_clusters, minimap2, trim_from_paf,
│   │                               #   checkv, pick_best_trim, fetch_blast_input_seqs
│   └── blast_trim_tools.nf         # Step 4 processes: blastn, blastani, trim_genomes (blast mode),
│                                   #   merge_blast_trim, checkv, completeness_filter
├── subworkflows/
│   ├── cluster_trim.nf             # Step 3: minimap2-based cluster trimming subworkflow
│   └── blast_trim.nf               # Step 4: BLAST-mode reference trimming subworkflow
├── bin/                            # Python scripts (auto-added to PATH by Nextflow)
│   ├── parse_clusters.py           # Extract rank1/2/3 pairs from vclust clusters
│   ├── trim_from_paf.py            # Select best alignment per cluster from PAF, write BED
│   ├── pick_best_trim.py           # Split trimmed output by CheckV completeness
│   ├── trim_genomes.py             # Nucmer-based trimming for blast mode (step 4)
│   ├── blastani_nayfach.py         # Compute ANI from blastn tabular output
│   └── completeness_filter.py      # Filter sequences by CheckV completeness threshold
├── envs/
│   └── wvdb_build.yml              # Conda environment for all tools
└── test/
    └── data/                       # Small test inputs for local/cluster validation
```

---

## Dependencies

All tools are managed via the conda environment in `envs/wvdb_build.yml`.

| Tool | Version | Notes |
|---|---|---|
| Nextflow | ≥ 23.10 | Requires Java 11+ |
| seqkit | 2.9.0 | |
| vclust | 1.3.1 | |
| minimap2 | ≥ 2.26 | Cluster-mode trimming (step 3) |
| blastn | 2.16.0+ | Via `blast` conda package |
| checkv | latest (1.1.1+) | |
| nucmer | 4.0.1 | BLAST-mode trimming (step 4); via `mummer4` conda package |
| samtools | ≥ 1.18 | Provides `bgzip` for PAF compression |
| diamond | ≥ 2.0.9 | CheckV hard dependency |
| Python | 3.11 | biopython, pandas, numpy required by bin/ scripts |

---

## Installation

### 1. Create the conda environment

```bash
# Recommended: use mamba for faster dependency solving
conda install -n base -c conda-forge mamba

mamba env create -f envs/wvdb_build.yml
conda activate wvdb-build
```

### 2. Verify all tools resolve

Each tool uses a different version flag:

```bash
seqkit version                    # seqkit v2.9.0
vclust -v                         # vclust 1.3.1
minimap2 --version                # 2.26+
blastn -version                   # blastn: 2.16.0+
checkv -h 2>&1 | head -1          # checkv 1.1.x
nucmer --version                  # 4.0.1
nextflow -v                       # nextflow version 23.x
bgzip --version 2>&1 | head -1    # bgzip 1.18+
```

Verify Python helper scripts:

```bash
for script in parse_clusters.py trim_from_paf.py pick_best_trim.py \
              trim_genomes.py blastani_nayfach.py completeness_filter.py; do
    python bin/$script --help > /dev/null 2>&1 \
        && echo "OK: $script" \
        || echo "CHECK: $script — try running with no args"
done
```

### 3. Validate the pipeline DAG (no data needed)

```bash
cd wvdb-build
nextflow run main.nf -stub -profile local
```

---

## Running the pipeline

### SLURM cluster (production)

```bash
conda activate wvdb-build
cd wvdb-build

nextflow run main.nf -profile slurm,conda \
    --fastqdir /path/to/fastq \
    --outdir   /path/to/results \
    --checkvdb /path/to/checkv-db-v1.5 \
    --refseq_ev_blastdb /path/to/esviritu_plus_refseq_virus.fna \
    --imgvr_blastdb     /path/to/IMGVR5_UViG.fna \
    -resume
```

The `-profile slurm,conda` combination uses SLURM as the executor and
activates the wvdb-build conda environment for every submitted job.
This is required because SLURM jobs run in a fresh shell that does not
inherit the interactive `conda activate`.

The `-resume` flag restarts from the last successful checkpoint if any
process fails.

### Local development / macOS

BLAST-mode trimming (step 4) requires large reference databases that are
not practical to run locally. For local development, run steps 1–3 and 5
by temporarily commenting out `BLAST_TRIM` in `main.nf`.

```bash
conda activate wvdb-build
cd wvdb-build

nextflow run main.nf -profile local \
    --fastqdir /path/to/fastq \
    --outdir   /path/to/results \
    --checkvdb /path/to/checkv-db-v1.5 \
    --refseq_ev_blastdb /path/to/small_test_db.fna \
    --imgvr_blastdb     /path/to/small_test_db.fna \
    --max_memory '32 GB' \
    -resume
```

### Test profile

Place test FASTAs under `test/data/` following the expected layout:

```
test/data/
└── <sample>/
    └── assembly/
        └── <run>/
            ├── filtered_virus.fasta
            └── filtered_provirus.fasta
```

```bash
nextflow run main.nf -profile test,conda \
    --checkvdb /path/to/checkv-db-v1.5 \
    --refseq_ev_blastdb /path/to/esviritu_plus_refseq_virus.fna \
    --imgvr_blastdb     /path/to/IMGVR5_UViG.fna
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--fastqdir` | required | Top-level input dir; expects `<sample>/assembly/<run>/filtered_*.fasta` |
| `--outdir` | required | Output directory |
| `--checkvdb` | required | CheckV database directory (e.g. `checkv-db-v1.5`) |
| `--refseq_ev_blastdb` | required | RefSeq virus + EsViritu BLAST db (`.fna` + all index files in same dir) |
| `--imgvr_blastdb` | required | IMG-VR v5 (MetaVR) BLAST db (`.fna` + all index files in same dir) |
| `--ani` | 0.95 | ANI threshold for vclust clustering |
| `--qcov` | 0.85 | Query coverage threshold for vclust |
| `--completeness` | 90 | CheckV completeness threshold (%) |
| `--threads` | 100 | Thread count for high-CPU processes |
| `--max_memory` | `'64 GB'` | Memory cap; reduce for local runs (e.g. `'32 GB'`) |
| `--conda_env` | auto | Path to conda env yml or prefix; overrides default in `conda` profile |

> **BLAST database note:** BLAST index files must reside in the same directory
> as the `.fna` file. The pipeline passes database paths as strings (not staged
> files) to preserve index file access. Do not symlink individual files.

---

## Output structure

```
outdir/
├── 1_collect/
│   ├── filtered_all.fasta              # all input genomes merged
│   ├── filtered_all.length.txt
│   ├── filtered_virus_all.fasta
│   └── filtered_provirus_all.fasta
├── 2_clustered/
│   ├── vclust_clusters.tsv             # cluster assignments (object<TAB>cluster)
│   └── vclust_ani.ids.tsv              # sequence ID/length index (produced by vclust align)
├── 3_cluster_trim/
│   ├── pairs12.tsv                     # rank1 vs rank2 pairs
│   ├── pairs13.tsv                     # rank1 vs rank3 pairs
│   ├── pairs23.tsv                     # rank2 vs rank3 pairs
│   ├── singletons.txt                  # 1-member clusters → step 4
│   ├── candidate_ids.txt               # rank1/2/3 IDs for minimap2
│   ├── trimming_candidates.fasta       # extracted rank1/2/3 sequences
│   ├── alignments.paf.gz               # full minimap2 PAF (bgzipped) for network analysis
│   ├── trimmed.fasta                   # best-trimmed sequence per cluster
│   ├── trimming.bed                    # BED coordinates used for trimming
│   ├── checkv/                         # CheckV quality assessment
│   ├── complete_reps.fasta             # complete trimmed reps → step 5
│   └── blast_trim_ids.txt              # incomplete cluster IDs → step 4
├── 4_blast_trim/
│   ├── blast_trim_input.fasta          # singletons + incomplete from step 3
│   ├── refseq_ev/                      # blastn + blastani + trim outputs (refseq_ev db)
│   ├── metavr/                         # blastn + blastani + trim outputs (metavr db)
│   ├── blast_trim_merged.fasta         # merged trimmed output
│   ├── checkv/
│   └── cluster_reps_complete.fasta     # complete reps → step 5
└── 5_reclustered/
    ├── recluster_input.fasta           # merged complete reps from steps 3 + 4
    ├── vclust_clusters.tsv
    ├── vclust_centroids.txt
    └── vclust_centroids.fasta          # ← FINAL OUTPUT: non-redundant genome set
```

---

## Pipeline overview

```
Input FASTAs (filtered_virus.fasta + filtered_provirus.fasta per sample)
    │
    ▼
Step 1: Collect & merge (COLLECT_GENOMES)
    │
    ▼
Step 2: Cluster (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER)
    │
    ▼
Step 3: Cluster trimming (CLUSTER_TRIM subworkflow)
    │   parse_clusters.py → rank1/2/3 pairs + candidate sequences
    │   minimap2 all-vs-all on candidates → alignments.paf
    │   trim_from_paf.py → pick best alignment per cluster (pairs12/13/23)
    │                       by alignment block length → BED → trimmed.fasta
    │   CheckV → pick_best_trim.py
    │       complete → step 5
    │       incomplete + singletons → step 4
    │
    ▼
Step 4: BLAST-mode trimming (BLAST_TRIM subworkflow)
    │   singletons + incomplete from step 3
    │   blastn × 2 (refseq_ev + metavr, parallel) → blastani → trim_genomes.py
    │   merge (refseq_ev results + metavr-only) → CheckV → completeness_filter
    │       complete → step 5
    │
    ▼
Step 5: Recluster (VCLUST_PREFILTER → VCLUST_ALIGN → VCLUST_CLUSTER → GET_CENTROIDS)
    │
    ▼
vclust_centroids.fasta — final non-redundant viral genome set
```

---

## Development notes

### Running from the correct directory

Always run Nextflow from inside `wvdb-build/`:

```bash
cd /path/to/repo/wvdb-build
nextflow run main.nf ...
```

Nextflow resolves `./modules/`, `./subworkflows/`, and `./bin/` relative to
`main.nf`. Running from the repo root will produce "invalid include source" errors.

### Conda environment on the cluster

SLURM jobs run in a fresh shell and do not inherit `conda activate` from the
login shell. Always use `-profile slurm,conda` (not just `-profile slurm`) for
cluster runs. The `conda` profile sets `process.conda` so every submitted job
activates the environment before running.

### BLAST database staging

BLAST databases are passed as `val` strings (not `path` inputs) to prevent
Nextflow from staging only the `.fna` file and leaving the index files behind.
All index files must be present in the same directory as the `.fna` file.

### Updating the conda environment

```bash
# Edit envs/wvdb_build.yml, then:
mamba env update -n wvdb-build -f envs/wvdb_build.yml --prune
nextflow run main.nf -stub -profile local   # confirm DAG still validates
git add envs/wvdb_build.yml && git commit -m "feat: update conda env"
```

### PAF output for network analysis

`3_cluster_trim/alignments.paf.gz` contains the full minimap2 all-vs-all
alignment of rank1/2/3 candidate sequences. This file can be used to build
a pairwise similarity network of cluster representatives for downstream
visualization (e.g. Cytoscape, gephi) or to re-examine clustering thresholds
without rerunning minimap2.

### task.ext.publish_dir pattern

Processes in `cluster_trim_tools.nf` and `blast_trim_tools.nf` use
`publishDir { "${params.outdir}/${task.ext.publish_dir}" }` with a closure.
The closure form is required because `task.ext` is only populated at runtime.
Values are set via `withName` selectors in `nextflow.config`.
