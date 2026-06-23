# wf3-pipeline

Viral genome clustering and reference-guided trimming pipeline, converted from `run_wf3_v2.sh`.

## Dependencies

| Tool | Version |
|---|---|
| Nextflow | ≥ 23.10 |
| seqkit | 2.9.0 |
| vclust | 1.3.1 |
| blastn | 2.16.0+ |
| checkv | 1.0.3 |
| nucmer | 4.0.1 |

Python helper scripts (place in `bin/` or ensure they are on `$PATH`):
- `blastani_nayfach.py`
- `parse_clusters_v2.py`
- `trim_genomes.py`
- `completeness_filter.py`

## Quick start

```bash
# 1. Install Nextflow (needs Java 11+)
curl -s https://get.nextflow.io | bash

# 2. Validate DAG without executing (no tools needed)
nextflow run main.nf -stub

# 3. Run locally
nextflow run main.nf \
    --fastqdir /path/to/fastq \
    --outdir   /path/to/results \
    --checkvdb /path/to/checkv-db-v1.5 \
    --refseq_ev_blastdb /path/to/esviritu_plus_refseq_virus.fna \
    --imgvr_blastdb     /path/to/IMGVR5_UViG.fna

# 4. Run on SLURM
nextflow run main.nf -profile slurm \
    --fastqdir /path/to/fastq \
    --outdir   /path/to/results \
    ...
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--fastqdir` | required | Top-level wf2 output dir (`<sample>/assembly/<run>/filtered_*.fasta`) |
| `--outdir` | required | Output directory |
| `--checkvdb` | required | CheckV database directory |
| `--refseq_ev_blastdb` | required | RefSeq + EsViritu combined BLAST db |
| `--imgvr_blastdb` | required | IMG-VR v5 BLAST db |
| `--ani` | 0.95 | ANI threshold for clustering |
| `--qcov` | 0.85 | Query coverage threshold |
| `--completeness` | 90 | CheckV completeness threshold (%) |
| `--threads` | 100 | Threads (used as default cpus for high-cpu processes) |

## Output structure

```
outdir/
├── 1_filtered/
│   ├── filtered_all.fasta
│   ├── filtered_all.length.txt
│   ├── filtered_virus_all.fasta
│   └── filtered_provirus_all.fasta
├── 2_clustered/
│   ├── vclust_clusters.tsv
│   └── vclust_ani.ids.tsv
├── 3_branchA/          (trimmed, checkv, complete/incomplete)
├── 4_branchB/
├── 5_branchC/          (TODO)
└── 6_reclustered/
    ├── vclust_clusters.tsv
    ├── vclust_centroids.txt
    └── vclust_centroids.fasta  ← final output
```

## Status

- [x] Step 1 — Collect genomes
- [x] Step 2 — vclust clustering
- [x] Steps 3/4 — Branch A and B (shared TRIM_AND_FILTER subworkflow)
- [ ] Step 5 — Branch C (BLAST-mode trimming) — `subworkflows/branchC.nf` TODO
- [x] Step 6 — Reclustering
