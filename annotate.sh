
###run checkV###
workdir=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3
checkv=/usr/WS1/mlbiomon/software/mamba/envs/checkv/bin/checkv
checkvdb=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/checkv_db/checkv-db-v1.5
$checkv end_to_end -d $checkvdb -t 100 $workdir/votus_final.fasta $workdir/checkv

###run genoNomad###
genomad_db=/p/vast1/mlbiomon/ref_data/genomad_db_1.9/genomad_db
in=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/votus_final.fasta
out=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/genomad/
genomad end-to-end --cleanup $in $out $genomad_db -t 120

###run RdRPCATCH###
rdrpcatch=/usr/WS1/mlbiomon/software/rdrpcatch/bin/rdrpcatch
db=/p/vast1/mlbiomon/ref_data/rdrp_catch_db
fasta=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/votus_final.fasta
outdir=/p/vast1/mlbiomon/analysis/securebio/viral_genomes_analysis/clusters_v3/rdrpcatch/rdrpcatch_out

$rdrpcatch scan -i $fasta -o $outdir -db-dir $db --cpus 128 -seq_type nuc


###run BLASTn###
# should allow user to specify desired databases, ideally in a config because command line could get long
# user should specify a name for each database so that hits can be reported clearly
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

###run BLASTp or DIAMOND###
#use predicted proteins from prodigal-gv output by geNomad
#run against NCBI-nr

###run RNAVirHost###
# must run after geNomad and checkV so that order.csv can be generated as input first
bin/run_rnavirhost.py # see options based on this script
# not fully tested

###Merge annotations###
bin/merge_annotations.py # see options based on this script
# this script initially had run_rnavirhost.py embedded within, may still have leftovers
# not fully tested

###Guess host by ensemble approach###
bin/guess_host.py # see options based on this script
# this is the roughest script, currently requires LLM access but this should be optional
# I copy-pasted to/from LLM but the prompt calls should be tied in directly (with API key)