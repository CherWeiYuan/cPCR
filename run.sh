# Run primer design

dir=~/Desktop/pcr/
csv_file=$dir/trial/input/target_exons.csv
preset=sPCR_gencode
prefix=trial
outdir=$dir/trial
threads=16


mamba activate pcr
pcr \
--target_exons_csv $csv_file \
--preset $preset \
--prefix $prefix \
--outdir $outdir \
--genome_fasta $dir/genome/GRCh38.p14.genome.fasta \
--gtf_file $dir/genome/gencode.v47.annotation.gtf \
--off_target_fasta $dir/transcriptome/gencode.v47.transcripts.fasta \
--blast_db_name gencode.v47.transcripts \
--blastdb_dir $dir/blastdb \
--constitutive_threshold 0.8 \
--polymerase dreamtaq \
--delete_xml \
--threads $threads
mamba deactivate