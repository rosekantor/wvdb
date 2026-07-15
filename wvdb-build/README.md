# wvdb_build

Viral genome database build pipeline. Collects viral and proviral genome
assemblies across samples, clusters them by sequence identity, trims chimeric
assemblies using nucmer pairwise alignments, assesses completeness with CheckV,
and produces a non-redundant set of representative genomes.

This pipeline lives in the `wvdb-build/` subdirectory of the wvdb repository,
alongside `wvdb-annotate/` and other workflows.

---

## Pipeline overview

```mermaid
flowchart TD
    A[/"Input FASTAs\nfiltered_virus + filtered_provirus\nper sample"/]
    A --> B

    B["Step 1: Collect\nCOLLECT_GENOMES\nfiltered_all.fasta"]
    B --> C

    C["Step 2: Cluster\nVCLUST_PREFILTER → ALIGN → CLUSTER\nvclust_clusters.tsv"]
    C -->|singletons| D
    C -->|multi-member clusters| E

    subgraph CT ["Step 3: Cluster trim"]
        E["PARSE_CLUSTERS + EXTRACT_TRIMMING_SEQS\nrank1/2/3 pairs + candidates.fasta"]
        E --> F12["TRIM_GENOMES trim12\nrank1 vs rank2"]
        E --> F13["TRIM_GENOMES trim13\nrank1 vs rank3"]
        E --> F23["TRIM_GENOMES trim23\nrank2 vs rank3"]
        F12 --> G["MERGE_TRIMMED\nsuffixed IDs → all_trimmed.fasta"]
        F13 --> G
        F23 --> G
        G --> H["CHECKV (one run)\nquality_summary.tsv"]
        H --> I["PICK_BEST_TRIM\nbest per rank1 by completeness + length"]
    end

    I -->|"complete ✓"| R
    I -->|incomplete| D

    subgraph BT ["Step 4: BLAST trim (optional)"]
        D["COLLECT_BLAST_INPUT\nsingletons + cluster-trim incomplete"]
        D --> J["BLASTN initial db"]
        D --> K["BLASTN secondary db\n(if run_secondary_blast=true)"]
        J --> L["SELECT_BEST_BLAST_HIT\nprefer initial; secondary if no initial hit"]
        K --> L
        L -->|"no hit in either db"| U
        L -->|"has hit"| M["TRIM_GENOMES_BLAST + CHECKV\nnucmer trim → completeness filter"]
    end

    M -->|"complete ✓"| R
    M -->|incomplete| U

    U[("COLLECT_UNVALIDATED\nunvalidated_genomes.fasta\nunvalidated_report.tsv")]

    R["Step 5: Recluster\nVCLUST complete reps from steps 3+4\n→ vclust_centroids.fasta"]
    R --> S

    S[/"vclust_centroids.fasta\nfinal non-redundant viral genome set"/]
    S --> T["PIPELINE_SUMMARY\npipeline_summary.md + .tsv"]

    style CT fill:none,stroke:#888,stroke-dasharray:5 3
    style BT fill:none,stroke:#888,stroke-dasharray:5 3
    style U fill:#faeeda,stroke:#ba7517
    style S fill:#e1f5ee,stroke:#0f6e56
```

### Sequence routing summary

| Sequence type | After step 3 | After step 4 | Final destination |
|---|---|---|---|
| Multi-member cluster, trimmed complete | → recluster | — | `vclust_centroids.fasta` |
| Multi-member cluster, trimmed incomplete | → blast_trim | complete → recluster | `vclust_centroids.fasta` |
| Multi-member cluster, trimmed incomplete | → blast_trim | incomplete → unvalidated | `unvalidated_genomes.fasta` |
| Singleton | → blast_trim | complete → recluster | `vclust_centroids.fasta` |
| Singleton | → blast_trim | no hit or incomplete | `unvalidated_genomes.fasta` |
| Any | — | — (run_initial_blast=false) | `unvalidated_genomes.fasta` |

### Step 3 trimming detail

All three trim passes (trim12, trim13, trim23) run in parallel using the
pre-indexed `trimming_candidates.fasta`. Their outputs are merged with
suffixed IDs (`<rank1_id>_trim12` etc.) so a **single CheckV run** covers
all three. `PICK_BEST_TRIM` then selects the best result per rank1 sequence
by completeness (primary) and length (tiebreaker), with trim12 preferred
on a complete tie.

---

## Repository structure

```
wvdb-build/
├── main.nf                         # Entry point; wires all steps
├── nextflow.config                 # Params, profiles (local, slurm, slurm_nomem, cluster, conda, test)
├── conf/
│   └── base.config                 # Per-process CPU/memory resource labels
├── modules/
│   ├── collect.nf                  # Step 1: collect and merge input FASTAs
│   ├── vclust.nf                   # Steps 2 & 5: vclust prefilter/align/cluster + centroids
│   ├── cluster_trim_tools.nf       # Step 3 processes: parse_clusters, extract_seqs,
│   │                               #   trim_genomes (×3), merge_trimmed, checkv, pick_best_trim,
│   │                               #   fetch_blast_input_seqs, completeness_filter
│   ├── blast_trim_tools.nf         # Step 4 processes: blastn, blastani, select_best_hit,
│   │                               #   trim_genomes_blast, merge_trimmed, checkv,
│   │                               #   completeness_filter, collect_unvalidated
│   └── summary.nf                  # PIPELINE_SUMMARY: sequence counts report
├── subworkflows/
│   ├── cluster_trim.nf             # Step 3: parallel nucmer trimming → single CheckV
│   └── blast_trim.nf               # Step 4: BLAST-mode reference trimming (optional)
├── bin/                            # Python scripts (auto-added to PATH by Nextflow)
│   ├── parse_clusters.py           # Extract rank1/2/3 pairs from vclust clusters
│   ├── trim_genomes.py             # Nucmer-based trimming (cluster mode + BLAST mode)
│   ├── merge_trimmed.py            # Combine trim12/13/23 FASTAs with suffixed IDs for CheckV
│   ├── pick_best_trim.py           # Select best trimming per rank1 by completeness + length
│   ├── select_best_blast_hit.py    # Route queries to initial or secondary db
│   ├── blastani_nayfach.py         # Compute pairwise ANI from blastn tabular output
│   ├── completeness_filter.py      # Filter sequences by CheckV completeness threshold
│   ├── collect_unvalidated.py      # Gather unvalidated seqs with reason report
│   └── pipeline_summary.py         # Sequence counts at each step → TSV + markdown report
├── envs/
│   └── wvdb_build.yml              # Conda environment for all tools
└── test/
    └── data/                       # Small test inputs (sample/assembly/run/ structure)
```

---

## Dependencies

All tools are managed via the conda environment in `envs/wvdb_build.yml`.

| Tool | Version | Notes |
|---|---|---|
| Nextflow | ≥ 23.10 | Requires Java 11+ |
| seqkit | 2.9.0 | |
| vclust | 1.3.1 | |
| blastn | 2.16.0+ | Via `blast` conda package |
| checkv | ≥ 1.0.3 | |
| nucmer | 4.0.1 | All trimming steps; via `mummer4` conda package |
| diamond | ≥ 2.0.9 | CheckV hard dependency |
| Python | 3.11 | biopython, pandas, numpy required by bin/ scripts |

---

## Installation

### 1. Create the conda environment

```bash
# Recommended: use mamba for faster solving
conda install -n base -c conda-forge mamba
mamba env create -f envs/wvdb_build.yml

# If mamba fails to solve, fall back to conda:
# conda env create -f envs/wvdb_build.yml

conda activate wvdb-build
```

### 2. Verify all tools resolve

```bash
seqkit version                    # seqkit v2.9.0
vclust -v                         # vclust 1.3.1
blastn -version                   # blastn: 2.16.0+
checkv -h 2>&1 | head -1          # checkv 1.x
nucmer --version                  # 4.0.1
nextflow -v                       # nextflow version 23.x
```

Verify Python helper scripts:

```bash
for script in parse_clusters.py trim_genomes.py merge_trimmed.py pick_best_trim.py \
              select_best_blast_hit.py blastani_nayfach.py completeness_filter.py \
              collect_unvalidated.py pipeline_summary.py; do
    python bin/$script --help > /dev/null 2>&1 \
        && echo "OK: $script" \
        || echo "CHECK: $script"
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
    --fastqdir          /path/to/fastq \
    --outdir            /path/to/results \
    --checkvdb          /path/to/checkv-db-v1.5 \
    --initial_blastdb   /path/to/esviritu_plus_refseq_virus.fna \
    --secondary_blastdb /path/to/IMGVR5_UViG.fna \
    --run_secondary_blast true \
    -resume
```

The `-profile slurm,conda` combination uses SLURM as the executor and
activates the conda environment for every submitted job. This is required
because SLURM jobs do not inherit the interactive `conda activate`.

For clusters where `DefMemPerNode=UNLIMITED` (memory not tracked per job),
use the `slurm_nomem` profile instead of `slurm`. This is a custom profile
defined in `nextflow.config` that omits the `--mem` flag from sbatch
submissions, letting SLURM manage memory allocation automatically:

```bash
nextflow run main.nf -profile slurm_nomem,conda ...
```

For the dedicated 128-CPU / 2TB cluster nodes, use the `cluster` profile
which requests full exclusive node access and sets threads/memory
accordingly:

```bash
nextflow run main.nf -profile cluster,conda \
    --fastqdir          /path/to/fastq \
    --outdir            /path/to/results \
    --checkvdb          /path/to/checkv-db \
    --initial_blastdb   /path/to/esviritu_plus_refseq_virus.fna \
    -resume
```

### Initial BLAST database only (no secondary)

```bash
nextflow run main.nf -profile slurm,conda \
    --fastqdir        /path/to/fastq \
    --outdir          /path/to/results \
    --checkvdb        /path/to/checkv-db \
    --initial_blastdb /path/to/esviritu_plus_refseq_virus.fna \
    -resume
```

### Skipping BLAST trim entirely

When `--run_initial_blast false`, all singletons and cluster-trim
incomplete sequences go directly to `unvalidated/` with no BLAST search:

```bash
nextflow run main.nf -profile slurm,conda \
    --fastqdir          /path/to/fastq \
    --outdir            /path/to/results \
    --checkvdb          /path/to/checkv-db \
    --run_initial_blast false \
    -resume
```

### Local development / macOS

```bash
conda activate wvdb-build
cd wvdb-build

nextflow run main.nf -profile local \
    --fastqdir          /path/to/fastq \
    --outdir            /path/to/results \
    --checkvdb          /path/to/checkv-db \
    --run_initial_blast false \
    --max_memory        '32 GB' \
    -resume
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--fastqdir` | required | Top-level input dir; expects `<sample>/assembly/<run>/filtered_*.fasta` |
| `--outdir` | required | Output directory |
| `--checkvdb` | required | CheckV database directory (e.g. `checkv-db-v1.5`) |
| `--initial_blastdb` | null | Primary BLAST reference db — required if `run_initial_blast=true` |
| `--secondary_blastdb` | null | Secondary BLAST reference db — required if `run_secondary_blast=true` |
| `--run_initial_blast` | true | Search `initial_blastdb` for singletons + incomplete sequences |
| `--run_secondary_blast` | false | Also search `secondary_blastdb`; initial db hit preferred when both match |
| `--ani` | 0.95 | ANI threshold for vclust clustering |
| `--qcov` | 0.85 | Query coverage threshold for vclust |
| `--completeness` | 90 | CheckV completeness threshold (%) |
| `--threads` | 100 | Thread count for high-CPU processes |
| `--max_memory` | `'64 GB'` | Memory cap for high/medium-CPU processes |
| `--max_memory_low` | `'16 GB'` | Memory cap for low-CPU processes |
| `--max_time` | `'24 h'` | Maximum runtime for any process |
| `--conda_env` | auto | Conda env yml path or prefix; overrides default in `conda` profile |

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
│   └── vclust_ani.ids.tsv              # sequence ID/length file from vclust align
├── 3_cluster_trim/
│   ├── pairs12.tsv                     # rank1 vs rank2 pairs
│   ├── pairs13.tsv                     # rank1 vs rank3 pairs
│   ├── pairs23.tsv                     # rank2 vs rank3 pairs
│   ├── singletons.txt                  # 1-member clusters → step 4
│   ├── candidate_ids.txt               # rank1/2/3 IDs for nucmer
│   ├── trimming_candidates.fasta       # extracted rank1/2/3 sequences
│   ├── trim12/                         # trim12 nucmer outputs
│   ├── trim13/                         # trim13 nucmer outputs
│   ├── trim23/                         # trim23 nucmer outputs
│   ├── all_trimmed.fasta               # merged with suffixed IDs for CheckV
│   ├── checkv/                         # CheckV quality assessment (single run)
│   ├── complete_reps.fasta             # complete trimmed reps → step 5
│   └── blast_trim_ids.txt              # incomplete cluster IDs → step 4
├── 4_blast_trim/                       # only created if run_initial_blast=true
│   ├── blast_trim_input.fasta          # singletons + incomplete from step 3
│   ├── initial/                        # blastn + blastani + trim outputs (initial db)
│   ├── secondary/                      # blastn + blastani + trim outputs (secondary db)
│   ├── blast_trimmed.fasta             # merged trimmed output
│   ├── checkv/                         # CheckV quality assessment
│   └── cluster_reps_complete.fasta     # complete reps → step 5
├── unvalidated/
│   ├── unvalidated_genomes.fasta       # untrimmed seqs: no BLAST hit or incomplete after trim
│   └── unvalidated_report.tsv          # per-seq reason: no_blast_hit | incomplete_after_trimming
├── 5_reclustered/
│   ├── recluster_input.fasta           # merged complete reps from steps 3 + 4
│   ├── vclust_clusters.tsv
│   ├── vclust_centroids.txt
│   └── vclust_centroids.fasta          # ← FINAL OUTPUT: non-redundant genome set
├── pipeline_summary.tsv                # sequence counts at each step (machine-readable)
└── pipeline_summary.md                 # sequence counts report (human-readable)
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
login session. Always use `-profile slurm,conda` for cluster runs.

### BLAST database staging

BLAST databases are passed as `val` strings (not `path` inputs) to prevent
Nextflow from staging only the `.fna` file away from its index files. All
index files must remain in the same directory as the database file.

### Memory configuration across clusters

| Cluster config | Profile to use | Notes |
|---|---|---|
| `DefMemPerCPU` set | `slurm` | Explicit `--mem` requests work normally |
| `DefMemPerNode=UNLIMITED` | `slurm_nomem` | Custom profile — omits `--mem` flag |
| 128-CPU / 2TB exclusive nodes | `cluster` | Custom profile — `--exclusive`, 128 CPUs, memory=null |

All three are custom profiles defined in `nextflow.config` under the `profiles {}`
block — not built-in Nextflow terms. Override memory caps at runtime:

```bash
nextflow run main.nf --max_memory '120 GB' --max_memory_low '8 GB' ...
```

### Updating the conda environment

```bash
mamba env update -n wvdb-build -f envs/wvdb_build.yml --prune
# If mamba fails to solve:
# conda env remove -n wvdb-build
# conda env create -f envs/wvdb_build.yml
nextflow run main.nf -stub -profile local
git add envs/wvdb_build.yml && git commit -m "feat: update conda env"
```

### task.ext.publish_dir pattern

Processes in `cluster_trim_tools.nf` and `blast_trim_tools.nf` use
`publishDir { "${params.outdir}/${task.ext.publish_dir}" }` with a closure.
The closure form is required because `task.ext` is only populated at runtime.
Values are set via `withName` selectors in `nextflow.config`.
