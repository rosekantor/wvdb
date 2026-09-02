# wvdb_annotate

Viral genome annotation pipeline. Takes a FASTA file of viral genomes
(e.g. the `vclust_centroids.fasta` output of `wvdb_build`) and produces
completeness assessment, classification, RdRP detection, nucleotide and
protein similarity searches, host prediction, and a merged per-vOTU
annotation table.

This pipeline lives in the `wvdb-annotate/` subdirectory of the wvdb
repository, alongside `wvdb-build/`.

---

## Pipeline overview

```mermaid
flowchart TD
    A[/"--input_fasta\nvOTU FASTA"/]

    A --> B["CheckV\ncompleteness + quality"]
    A --> C["geNomad\nclassification + gene prediction"]
    A --> D["RdRPCATCH\nRdRP detection\n(run_rdrpcatch=true)"]
    A --> E["BLASTn × N databases\nnucleotide similarity\n(run_blastn=true)"]
    E --> F["BLASTani\npairwise ANI"]
    C -->|proteins.faa| P["CHARACTERIZE_PROTEINS\nDIAMOND × N + hmmsearch × N\n(run_diamond / run_hmmsearch)"]
    C --> H["RNAVirHost\nhost prediction\n(run_rnavirhost=true)"]
    D --> H

    B --> I["MERGE_ANNOTATIONS\nper-vOTU annotation table"]
    C --> I
    D --> I
    F --> I
    P --> I
    H --> I

    I --> J["CATEGORIZE_HOST_TERMS\nLLM categorization (cacheable)\n(run_guess_host=false)"]
    J --> L["GUESS_HOST\ndecision tree, no LLM dependency"]
    I --> K["ANNOTATION_SUMMARY\nannotation_summary.md + .tsv"]
    L --> K

    style A fill:#e1f5ee,stroke:#0f6e56
    style K fill:#e1f5ee,stroke:#0f6e56
```

All steps except RNAVirHost run in parallel from the input FASTA.
RNAVirHost runs sequentially after geNomad and RdRPCATCH complete, as it
requires their outputs to generate its consensus order classification
(the input FASTA itself is the canonical source of contig IDs — CheckV
is not part of this dependency).

### Step control

| Step | Default | Param to disable/enable |
|---|---|---|
| CheckV | always on | — |
| geNomad | always on | — |
| RdRPCATCH | on | `--run_rdrpcatch false` |
| BLASTn | on | `--run_blastn false` |
| DIAMOND | off | `--run_diamond true` |
| hmmsearch | off | `--run_hmmsearch true` |
| RNAVirHost | on | `--run_rnavirhost false` |
| GUESS_HOST | off | `--run_guess_host true` (requires `--llm_env_file`, `--llm_provider`, `--llm_model`) |

---

## Repository structure

```
wvdb-annotate/
├── main.nf                         # Entry point; wires all annotation steps
├── nextflow.config                 # Params, profiles (local, slurm, slurm_nomem, cluster, conda, test)
├── conf/
│   └── base.config                 # Per-process CPU/memory resource labels
├── modules/
│   ├── checkv.nf                   # CheckV quality assessment
│   ├── genomad.nf                  # geNomad classification + gene prediction
│   ├── rdrpcatch.nf                # RdRPCATCH RdRP detection (separate conda env)
│   ├── blastn_tools.nf             # BLASTN + BLASTANI (one process per database, parallel)
│   ├── protein_tools.nf            # DIAMOND, HMMSEARCH, PROTEIN_SUMMARY
│   ├── rnavirhost.nf               # RNAVirHost host prediction (separate conda env)
│   └── summary.nf                  # MERGE_ANNOTATIONS, PREPARE_ICTV,
│                                   #   CATEGORIZE_HOST_TERMS, GUESS_HOST, ANNOTATION_SUMMARY
├── subworkflows/
│   └── characterize_proteins.nf   # DIAMOND × N + hmmsearch × N in parallel
├── bin/                            # Python scripts (auto-added to PATH by Nextflow)
│   ├── blastani_nayfach.py         # Compute pairwise ANI from blastn tabular output
│   ├── genomad_utils.py            # Shared: normalize geNomad provirus-suffixed contig IDs
│   ├── merge_annotations.py        # Combine all annotation outputs per vOTU
│   ├── run_rnavirhost.py           # RNAVirHost wrapper + consensus order generation
│   ├── llm_utils.py                # Shared: multi-vendor LLM dispatcher (anthropic/openai/
│   │                               #   openai_compatible), raw HTTP requests, no SDK deps
│   ├── categorize_host_terms.py    # LLM categorization of host/isolation-source text (cacheable)
│   ├── guess_host.py               # Ensemble host determination decision tree (no LLM dependency)
│   ├── prepare_ictv.py             # Process ICTV family HTML → ictv_families.tsv
│   ├── parse_hmmsearch.py          # Parse hmmsearch --domtblout to clean TSV
│   ├── summarize_protein_hits.py   # Per-contig protein hit counts
│   └── annotation_summary.py       # Per-step annotation coverage → TSV + markdown report
├── envs/
│   ├── wvdb_annotate.yml           # Main conda environment
│   ├── wvdb_rdrpcatch.yml          # Separate env for RdRPCATCH (Python 3.12)
│   └── wvdb_rnavirhost.yml         # Separate env for RNAVirHost (pinned ML deps)
├── ref_data/
│   └── ictv_families.tsv           # Bundled ICTV family table (MSL Jan 2026)
└── test/
    └── data/
        └── test_votus.fasta        # Small test input
```

---

## Dependencies

Three conda environments are required due to conflicting Python version and
dependency constraints between tools.

### Main environment (`wvdb-annotate`)

| Tool | Version | Notes |
|---|---|---|
| Nextflow | ≥ 23.10 | Requires Java 11+ |
| checkv | ≥ 1.0.3 | |
| genomad | ≥ 1.8.0 | Requires separate database download |
| blastn | 2.16.0+ | Via `blast` conda package |
| diamond | ≥ 2.0.9 | Protein search |
| hmmer | ≥ 3.3.2 | hmmsearch for HMM profile searches |
| Python | 3.11 | biopython (SeqIO, Entrez), pandas, numpy, requests, lxml, python-dotenv |

### RdRPCATCH environment (`rdrpcatch`)

Separate environment required — RdRPCATCH needs Python 3.12.

### RNAVirHost environment (`rnavirhost`)

Separate environment required — RNAVirHost requires pinned versions of
pandas (2.0.3), scikit-learn (1.1.3), and xgboost (1.7.4) that conflict
with other tools.

---

## Installation

### 1. Create the main conda environment

```bash
mamba env create -f envs/wvdb_annotate.yml
# If mamba fails: conda env create -f envs/wvdb_annotate.yml
conda activate wvdb-annotate
```

### 2. Create the RdRPCATCH environment

RdRPCATCH requires Python 3.12 and cannot share the main environment.

```bash
# Option A — install from bioconda directly (recommended)
conda create -n rdrpcatch -c bioconda rdrpcatch

# Option B — create from the provided yml
mamba env create -f envs/wvdb_rdrpcatch.yml
```

Download the RdRPCATCH databases:

```bash
conda activate rdrpcatch
rdrpcatch databases --destination-dir /path/to/rdrp_catch_db
```

The pipeline activates the `rdrpcatch` environment automatically for the
`RDRPCATCH` process via a per-process `conda` directive in `nextflow.config`.
Override the path if installed at a non-standard prefix:

```bash
nextflow run main.nf --rdrpcatch_conda_env /path/to/conda/envs/rdrpcatch ...
```

### 3. Create the RNAVirHost environment

RNAVirHost requires pinned ML dependency versions that conflict with the
main environment.

```bash
mamba env create -f envs/wvdb_rnavirhost.yml
conda activate rnavirhost
rnavirhost --help   # confirm installation
```

The pipeline activates this environment automatically for the `RNAVIRHOST`
process. Override the path if needed:

```bash
nextflow run main.nf --rnavirhost_conda_env /path/to/conda/envs/rnavirhost ...
```

### 4. Download the geNomad database

```bash
conda activate wvdb-annotate
genomad download-database /path/to/genomad_db/
```

Alternatively, download manually from Zenodo:
https://github.com/apcamargo/genomad

### 5. Download the CheckV database

```bash
checkv download_database /path/to/checkv-db/
```

### 6. Prepare the ICTV family reference file

A processed copy is bundled in `ref_data/ictv_families.tsv` (MSL January 2026)
and used by default — no setup required for most users. To update, see
[Updating reference data](#updating-reference-data).

### 7. Build BLAST databases

For each FASTA in your databases CSV, build the BLAST index if not already done:

```bash
makeblastdb -in /path/to/db.fna -dbtype nucl -out /path/to/db.fna
```

All index files must reside in the same directory as the `.fna` file.

### 8. Validate the pipeline DAG

```bash
cd wvdb-annotate
nextflow run main.nf -stub -profile local \
    --run_blastn false
```

### 9. Set up an LLM provider for host prediction (optional)

Only needed if you plan to use `--run_guess_host true`. Create a `.env`
file — suggested location `wvdb-annotate/.env`, but any absolute path is
accepted via `--llm_env_file`:

```bash
# Add this to your .gitignore if not already present
echo ".env" >> .gitignore

# Pick the line matching your --llm_provider
echo "ANTHROPIC_API_KEY=sk-ant-..." > .env
echo "OPENAI_API_KEY=sk-..." > .env
echo "OPENAI_COMPATIBLE_API_KEY=..." > .env   # for an internal enterprise gateway
```

See [Host prediction (GUESS_HOST)](#host-prediction-guess_host) below for
provider options and the full set of required parameters.

---

## Database configuration

All annotation databases (BLASTn, DIAMOND, hmmsearch profiles) are specified
in a single unified CSV file:

```csv
type,name,path
blastn,IMGVR,/data/ref_data/IMG-VR_2025-12-02/IMGVR5_UViG.fna
blastn,CHVD,/data/ref_data/CHVD_tisza2021/CHVD_virus_sequences_v1.1.fasta
blastn,UHGV,/data/ref_data/UHGV_nayfach2025/nayfach2025_votus_hq_plus.fna
blastn,VIRE,/data/ref_data/VIRE/all_vire.fna
blastn,core_nt,/data/core_nt/core_nt_Jul25_filt
diamond,nr,/data/ref_data/nr.dmnd
hmm,pfam,/data/ref_data/Pfam-A.hmm
```

- `type` — `blastn`, `diamond`, or `hmm`
- `name` — short label used for output filenames and merged annotation columns
- `path` — full path to the database file

Pass this file via `--databases /path/to/databases.csv`. Rows are only
processed when their type's corresponding flag is true (`run_blastn`,
`run_diamond`, `run_hmmsearch`). A row present in the CSV but with its
type's flag set to false is silently skipped. To add a new database,
add a row and rerun with `-resume` — only the new database will be searched.

**`core_nt` is special:** BLASTn hits against the database named `core_nt`
trigger an Entrez lookup for host organism and isolation source metadata.
Requires `--entrez_email` to be set.

---

## Running the pipeline

### SLURM cluster (production)

```bash
conda activate wvdb-annotate
cd wvdb-annotate

nextflow run main.nf -profile cluster,conda \
    --input_fasta   /path/to/votus_final.fasta \
    --outdir        /path/to/results \
    --checkvdb      /path/to/checkv-db-v1.5 \
    --genomad_db    /path/to/genomad_db \
    --rdrpcatch_db  /path/to/rdrp_catch_db \
    --databases     /path/to/databases.csv \
    --entrez_email  user@institution.edu \
    -resume
```

> **Avoiding unnecessary conda environment recreation:** by default, `-profile
> conda` points each environment-sensitive process at a yml file
> (`envs/wvdb_annotate.yml`, `envs/wvdb_rdrpcatch.yml`, `envs/wvdb_rnavirhost.yml`).
> When Nextflow sees a **yml path**, it always creates a fresh environment from
> it and caches that build under `work/conda/` — it does not search your conda
> installation for an existing environment with the same name. If you have
> already created `wvdb-annotate`, `rdrpcatch`, and `rnavirhost` environments
> manually, pass their **existing paths** instead of relying on the yml
> defaults, so Nextflow activates them directly with no build step:
>
> ```bash
> conda env list   # find exact paths
>
> nextflow run main.nf -profile cluster,conda \
>     --conda_env             /path/to/conda/envs/wvdb-annotate \
>     --rdrpcatch_conda_env   /path/to/conda/envs/rdrpcatch \
>     --rnavirhost_conda_env  /path/to/conda/envs/rnavirhost \
>     --input_fasta   /path/to/votus_final.fasta \
>     --outdir        /path/to/results \
>     --checkvdb      /path/to/checkv-db-v1.5 \
>     --genomad_db    /path/to/genomad_db \
>     --rdrpcatch_db  /path/to/rdrp_catch_db \
>     --databases     /path/to/databases.csv \
>     --entrez_email  user@institution.edu \
>     -resume
> ```

### With protein characterization enabled

```bash
nextflow run main.nf -profile cluster,conda \
    --input_fasta   /path/to/votus_final.fasta \
    --outdir        /path/to/results \
    --checkvdb      /path/to/checkv-db-v1.5 \
    --genomad_db    /path/to/genomad_db \
    --databases     /path/to/databases.csv \
    --entrez_email  user@institution.edu \
    --run_diamond   true \
    --run_hmmsearch true \
    -resume
```

### With host prediction enabled

```bash
nextflow run main.nf -profile cluster,conda \
    --input_fasta   /path/to/votus_final.fasta \
    --outdir        /path/to/results \
    --checkvdb      /path/to/checkv-db-v1.5 \
    --genomad_db    /path/to/genomad_db \
    --rdrpcatch_db  /path/to/rdrp_catch_db \
    --databases     /path/to/databases.csv \
    --entrez_email  user@institution.edu \
    --run_guess_host true \
    --llm_env_file  "$(pwd)/.env" \
    --llm_provider  anthropic \
    --llm_model     claude-sonnet-5 \
    -resume
```

Using an internal enterprise LLM gateway instead of a public provider:

```bash
nextflow run main.nf -profile cluster,conda \
    ... (same required params as above) \
    --run_guess_host true \
    --llm_env_file  "$(pwd)/.env" \
    --llm_provider  openai_compatible \
    --llm_model     gpt-4.1 \
    --llm_base_url  https://your-internal-gateway.example.org/v1/chat/completions \
    -resume
```

On a second run, reuse previously-categorized terms (and any hand
corrections you made to them) instead of re-querying the LLM:

```bash
nextflow run main.nf -profile cluster,conda \
    ... \
    --run_guess_host true \
    --llm_env_file  "$(pwd)/.env" \
    --llm_provider  anthropic \
    --llm_model     claude-sonnet-5 \
    --host_dict_cache      /path/to/results/annotations/host_categorization/host_dict.tsv \
    --isolation_dict_cache /path/to/results/annotations/host_categorization/isolation_dict.tsv \
    -resume
```

### Minimal run (CheckV + geNomad only)

```bash
nextflow run main.nf -profile cluster,conda \
    --input_fasta      /path/to/votus.fasta \
    --outdir           /path/to/results \
    --checkvdb         /path/to/checkv-db \
    --genomad_db       /path/to/genomad_db \
    --run_rdrpcatch    false \
    --run_blastn       false \
    --run_rnavirhost   false \
    -resume
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_fasta` | required | Input vOTU FASTA file |
| `--outdir` | required | Output directory |
| `--checkvdb` | required | CheckV database directory |
| `--genomad_db` | required | geNomad database directory |
| `--rdrpcatch_db` | null | RdRPCATCH database directory (required if `run_rdrpcatch=true`) |
| `--databases` | null | Unified databases CSV (required if any of `run_blastn`, `run_diamond`, `run_hmmsearch` is true) |
| `--ictv_fam` | bundled | Path to processed ICTV families TSV |
| `--entrez_email` | null | Email for NCBI Entrez (required for `core_nt` host lookup) |
| `--run_rdrpcatch` | true | Scan for RNA-dependent RNA polymerase |
| `--run_blastn` | true | BLASTn against `blastn` rows in databases CSV |
| `--run_diamond` | false | DIAMOND protein search against `diamond` rows in databases CSV |
| `--run_hmmsearch` | false | hmmsearch against `hmm` rows in databases CSV |
| `--run_rnavirhost` | true | Host prediction (requires geNomad + RdRPCATCH; FASTA is the canonical contig ID source, CheckV not required) |
| `--run_guess_host` | false | Ensemble host determination (ICTV + RNAVirHost + LLM-categorized BLAST metadata) |
| `--llm_env_file` | null | Path to `.env` file containing the API key — required if `run_guess_host=true` |
| `--llm_provider` | null | `anthropic` \| `openai` \| `openai_compatible` — required if `run_guess_host=true`, no default |
| `--llm_model` | null | Model name, e.g. `claude-sonnet-5`, `gpt-4o` — required if `run_guess_host=true`, no default |
| `--llm_base_url` | null | Required only if `llm_provider=openai_compatible` (e.g. an internal enterprise LLM gateway URL) |
| `--llm_api_key_env_var` | null | Override which `.env` variable holds the API key (default: `ANTHROPIC_API_KEY` / `OPENAI_API_KEY` / `OPENAI_COMPATIBLE_API_KEY` depending on provider) |
| `--host_dict_cache` | null | Path to a previously-generated `host_dict.tsv` to reuse (skips LLM calls for already-categorized terms) |
| `--isolation_dict_cache` | null | Path to a previously-generated `isolation_dict.tsv` to reuse |
| `--guess_host_other_threshold` | 100 | Host categories with fewer than this many vOTUs are consolidated into `"other"` |
| `--blastn_evalue` | `1e-3` | BLASTn e-value threshold |
| `--blastn_max_targets` | 10 | BLASTn max target sequences per query |
| `--diamond_evalue` | `1e-5` | DIAMOND e-value threshold |
| `--diamond_min_bitscore` | 50 | DIAMOND minimum bitscore |
| `--diamond_max_targets` | 10 | DIAMOND max target sequences per query |
| `--diamond_taxonmap` | false | Include staxids in DIAMOND output (requires `--taxonmap` at db build) |
| `--hmm_evalue` | `1e-5` | hmmsearch e-value fallback when `--cut_tc` is unavailable |
| `--threads` | 100 | Thread count for high-CPU processes |
| `--max_memory` | `'64 GB'` | Memory cap for high/medium-CPU processes |
| `--max_memory_low` | `'16 GB'` | Memory cap for low-CPU processes |
| `--max_time` | `'24 h'` | Maximum runtime for any process |
| `--rdrpcatch_conda_env` | `envs/wvdb_rdrpcatch.yml` | Path to rdrpcatch env or yml. **yml path → creates new env each time**; existing env path → activates directly, no rebuild |
| `--rnavirhost_conda_env` | `envs/wvdb_rnavirhost.yml` | Path to rnavirhost env or yml. Same yml-vs-path behavior as above |

---

## Output structure

```
outdir/
├── checkv/
│   └── checkv_out/
│       ├── quality_summary.tsv
│       ├── completeness.tsv
│       └── contamination.tsv
├── genomad/
│   └── genomad_out/
│       ├── <name>_summary/<name>_virus_summary.tsv
│       └── <name>_annotate/<name>_proteins.faa
├── rdrpcatch/                      # only if run_rdrpcatch=true
│   └── rdrpcatch_out/
├── blastn/                         # only if run_blastn=true
│   ├── IMGVR/
│   │   ├── IMGVR.blastn.tsv
│   │   └── IMGVR.blastn.ani.tsv
│   └── ...                         # one subdirectory per blastn row in databases.csv
├── proteins/                       # only if run_diamond or run_hmmsearch=true
│   ├── diamond/
│   │   └── <name>/
│   │       └── <name>.diamond.tsv
│   ├── hmmsearch/
│   │   └── <name>/
│   │       ├── <name>.hmmsearch.tsv
│   │       └── <name>.hmmsearch.domtbl
│   └── protein_summary.tsv         # per-contig protein hit counts
├── rnavirhost/                     # only if run_rnavirhost=true
│   ├── rnavirhost_out/
│   │   └── result.csv
│   └── rnavirhost_consensus_orders.csv
├── annotations/
│   ├── merged_annotations.tsv      # per-vOTU annotation table (all tools)
│   │                               #   includes genomad_provirus, genomad_provirus_start,
│   │                               #   genomad_provirus_end columns — see Development notes
│   │                               #   ⚠ if run_guess_host=true, this is an INTERMEDIATE —
│   │                               #   final_annotations.tsv below is the true final output
│   ├── host_categorization/        # only if run_guess_host=true
│   │   ├── host_dict.tsv           # LLM-categorized host terms (cacheable via --host_dict_cache)
│   │   ├── isolation_dict.tsv      # LLM-categorized isolation-source terms
│   │   └── rejected_terms.tsv      # terms that failed response validation, left
│   │                               #   uncategorized rather than guessed — retried next run
│   ├── final_annotations.tsv       # only if run_guess_host=true — merged_annotations.tsv's
│   │                               #   columns plus host_final and all decision-tree/review
│   │                               #   columns. THIS is the final output when host prediction
│   │                               #   is enabled; merged_annotations.tsv above is superseded.
│   └── manual_review.tsv           # only if run_guess_host=true — ambiguous order-level
│                                   #   cases (e.g. Martellivirales) not auto-corrected;
│                                   #   these rows still appear in final_annotations.tsv too
├── annotation_summary.tsv          # per-step annotation counts (machine-readable)
└── annotation_summary.md           # annotation counts report (human-readable)
```

---

## Updating reference data

### ICTV family table

A processed copy is bundled in `ref_data/ictv_families.tsv` (MSL January 2026).
To update when a new ICTV release is available:

**Step 1 — Save the HTML page manually:**
1. Open https://ictv.global/virus-properties in your browser
2. Set "Items per page" to "All" (bottom of page) to load all families
3. File → Save Page As → save as HTML

**Step 2 — Process and commit:**
```bash
nextflow run main.nf --prepare_ictv true \
    --ictv_raw /path/to/Virus_Properties___ICTV.html \
    --outdir   wvdb-annotate/ref_data
git add ref_data/ictv_families.tsv
git commit -m "ref: update ICTV family table to MSL <version>"
```

> **Note on new host mappings:** If ICTV adds new host combination strings
> not in the mapping table, `prepare_ictv.py` will warn about unmapped values
> and their `host_ICTV_simple` will be null. Add the new combination to the
> `HOST_REASSIGN` list in `bin/prepare_ictv.py` and rerun.
> Valid categories: `archaea`, `bacteria`, `fungi`, `invertebrates`, `plants`,
> `protists`, `vertebrates`, `non-vertebrates`, `incl-vertebrates`.

---

## Development notes

### Running from the correct directory

Always run Nextflow from inside `wvdb-annotate/`:

```bash
cd /path/to/repo/wvdb-annotate
nextflow run main.nf ...
```

### Three conda environments

| Environment | Python | Key constraint | Tools |
|---|---|---|---|
| `wvdb-annotate` | 3.11 | main env | checkv, genomad, blast, diamond, hmmer, all bin/ scripts |
| `rdrpcatch` | 3.12 | Python 3.12 required | rdrpcatch, mmseqs2 |
| `rnavirhost` | any | pandas=2.0.3, sklearn=1.1.3, xgboost=1.7.4 pinned | rnavirhost, prodigal |

The per-process conda env overrides for `RDRPCATCH` and `RNAVIRHOST` are set
in `nextflow.config` via `withName` selectors. Nextflow activates the correct
environment automatically for each process — no manual switching required.

### Memory configuration across clusters

| Cluster config | Profile | Notes |
|---|---|---|
| `DefMemPerCPU` set | `slurm` | Explicit `--mem` requests work normally |
| `DefMemPerNode=UNLIMITED` | `slurm_nomem` | Custom profile — omits `--mem` flag |
| 128-CPU / 2TB exclusive nodes | `cluster` | Custom profile — `--exclusive`, 128 CPUs, memory=null |

All three are custom profiles defined in `nextflow.config` — not built-in
Nextflow terms.

### BLAST database staging

BLAST databases are passed as `val` strings (not `path` inputs) to prevent
Nextflow from staging only the `.fna` file away from its index files.

### hmmsearch trusted cutoffs

`HMMSEARCH` attempts `--cut_tc` first (trusted cutoffs embedded in profiles
such as Pfam and VOG). If the profile has no TC thresholds, it falls back to
`--domE` with the value from `--hmm_evalue`. The fallback is logged to stderr.

### geNomad provirus-suffixed contig IDs

geNomad appends `|provirus_<start>_<end>` to a contig's ID in
`<prefix>_virus_summary.tsv` when it identifies an integrated provirus
region within that contig — one row per contig, never both a plain-ID
row and a suffixed row for the same contig. This suffixed ID does not
match the FASTA header, RdRPCATCH's output, or CheckV's `contig_id`,
which all use the plain ID. Left unhandled, this produced duplicate rows
in RNAVirHost's `orders.csv` (an outer join across differently-keyed
sources), which RNAVirHost rejects outright since it requires exactly
one row per FASTA sequence.

`bin/genomad_utils.py` is shared by `run_rnavirhost.py` and
`merge_annotations.py` and normalizes the ID (`load_genomad_virus_summary()`),
adding `genomad_provirus`, `genomad_provirus_start`, and
`genomad_provirus_end` columns rather than discarding the information.
This case matters biologically, not just technically: CheckV does not
always flag the same region as contamination, especially when the same
integrated element is reproducibly assembled across multiple independent
samples — that pattern is more consistent with a genuine phage/phage-plasmid
carrying real host genes than a chimeric assembly artifact. These
sequences are kept as-is and flagged in `merged_annotations.tsv` for
review, not trimmed or discarded automatically.

### Entrez robustness

`merge_annotations.py`'s Entrez lookup for `core_nt` hits reads each
batch's response into a string (`handle.read()`) before parsing with
`SeqIO.parse(io.StringIO(...), "genbank")`, rather than passing the live
`Entrez.efetch()` handle directly to `SeqIO.parse()`. The latter is
fragile across Biopython versions — `Bio.File.as_handle()` can fail to
recognize the handle type and crash with a confusing `TypeError` deep
inside Biopython on an otherwise-normal response. Each batch also retries
up to 3 times and, if it still fails, is skipped with a warning rather
than crashing the entire `MERGE_ANNOTATIONS` task — at ~65k sequences and
100 IDs per batch, one bad batch out of hundreds should not cost all the
others.

### Host prediction (GUESS_HOST)

> **For the full decision-tree rationale, diagrams, and the trust
> ordering between ICTV/BLAST/RNAVirHost, see [GUESS_HOST.md](GUESS_HOST.md).**

Host determination is split into two Nextflow processes so the slow,
API-cost-incurring part is cached independently of the cheap,
deterministic part:

**`CATEGORIZE_HOST_TERMS`** (LLM, cacheable) — extracts unique
`hosts_ntBlastHit`/`isolation_source_ntBlastHit` free-text values from
`merged_annotations.tsv` (populated via the `core_nt` Entrez lookup) and
categorizes each into a fixed set of categories using an LLM. Terms
already present in a supplied `--host_dict_cache`/`--isolation_dict_cache`
are not re-queried. Every LLM response is strictly validated: a term is
only accepted if its returned category is an exact match to one of the
allowed values; missing, blank, or invalid responses are left
**uncategorized** rather than defaulted to a fallback category, and are
retried automatically on the next run (see `rejected_terms.tsv`).

**`GUESS_HOST`** (decision tree, no LLM dependency) — combines ICTV
taxonomy, RNAVirHost's prediction, and the categorized BLAST metadata
into a final `host_final` call per vOTU, following the original
decision-tree logic (ICTV/RNAVirHost agreement checks, fecal-source
flagging, bacteriophage-class corrections for unclassified Orders,
always-bacteria Orders like Norzivirales/Timlovirales). Orders that are
*usually* one host type but not reliably enough to auto-correct (e.g.
Martellivirales is usually but not always plant-associated) are written
to `manual_review.tsv` instead of being silently corrected. This process
has no LLM dependency and is safe to rerun freely — e.g. after tuning the
`REVIEW_ORDERS`/`ALWAYS_BACTERIA_ORDERS` lists in `bin/guess_host.py` —
without re-querying the LLM.

**LLM providers:** `bin/llm_utils.py` supports `anthropic`, `openai`, and
`openai_compatible` via raw HTTP requests (no vendor SDK dependency).
`openai_compatible` accepts any endpoint implementing the same
chat-completions request/response shape as OpenAI — this covers internal
enterprise LLM gateways (a common pattern: front multiple vendor models
behind one OpenAI-shaped endpoint), Azure OpenAI, and locally-hosted
servers. Provider and model must always be specified explicitly — there
is no default model, since model availability changes over time and a
silently-used default could go stale unnoticed.

**`.env` handling:** the API key is loaded directly from the filesystem
by `categorize_host_terms.py`, never staged as a Nextflow `path` input —
this keeps it out of `work/` and any trace/log directory Nextflow
manages. Suggested location is `wvdb-annotate/.env` (gitignored), but any
absolute path works via `--llm_env_file`.

### Pending scripts

`annotation_summary.py` has been tested end-to-end on a full production
run and is confirmed working. `guess_host.py` and `categorize_host_terms.py`
are both implemented and tested (decision-tree logic validated against
synthetic cases covering every branch; LLM caching validated to make zero
API calls when all terms are pre-cached). Nothing is currently pending.
