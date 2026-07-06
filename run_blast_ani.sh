#!/bin/bash
#SBATCH -t 00-24:00:00
#SBATCH -N 1
#SBATCH -o stdout_blastn_%j
#SBATCH -e stderr_blastn_%j
#SBATCH --job-name blastn


###BLASTn###
workdir=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3
query=$workdir/votus_final.fasta
blastdb_IMGVR=/p/vast1/mlbiomon/ref_data/IMG-VR_2025-12-02/IMGVR5_UViG.fna
blastdb_corent=/p/vast1/kpath/blastdb/core_nt/core_nt_Jul25_filt
blastdb_CHVD=/p/vast1/mlbiomon/ref_data/CHVD_tisza2021/CHVD_virus_sequences_v1.1.fasta
blastdb_UHGV=/p/vast1/mlbiomon/ref_data/UHGV_nayfach2025/nayfach2025_votus_hq_plus.fna
blastdb_VIRE=/p/vast1/mlbiomon/ref_data/VIRE/all_vire.fna
out_imgvr=$workdir/blast_searches/votus_final-vs-IMGVRv5.tsv
out_corent=$workdir/blast_searches/votus_final-vs-corentJuly25.tsv
out_chvd=$workdir/blast_searches/votus_final-vs-CHVD.tsv
out_uhgv=$workdir/blast_searches/votus_final-vs-UHGV.tsv
out_vire=$workdir/blast_searches/votus_final-vs-VIRE.tsv

blastn -query $query -db $blastdb_IMGVR -evalue 1e-3 -max_target_seqs 10 -num_threads 100 -out $out_imgvr \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
blastn -query $query -db $blastdb_CHVD -evalue 1e-3 -max_target_seqs 10 -num_threads 100 -out $out_chvd \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
blastn -query $query -db $blastdb_UHGV -evalue 1e-3 -max_target_seqs 10 -num_threads 100 -out $out_uhgv \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
blastn -query $query -db $blastdb_VIRE -evalue 1e-3 -max_target_seqs 10 -num_threads 100 -out $out_vire \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
blastn -query $query -db $blastdb_corent -evalue 1e-3 -max_target_seqs 10 -num_threads 100 -out $out_corent \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

for i in $workdir/blast_searches/*tsv; do
  out=$(echo $i | sed 's/tsv/ani.tsv/g')
  /p/vast1/mlbiomon/analysis/securebio/scripts/blastani_nayfach.py -i $i -o $out
done
