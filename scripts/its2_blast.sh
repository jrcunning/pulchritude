#!/bin/bash
#SBATCH --job-name=its2_targ
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --account=putnamlab
#SBATCH -o its2.out
#SBATCH -e its2.error

#this script is designed to pull the designated ITS2 type from whole genome data to
#validate the bbsplit results 

source ~/.bashrc
conda activate repeatmasker 
module load parallel 

mkdir /data/putnamlab/tconn/apul_reseq/symits2
cd /data/putnamlab/tconn/apul_reseq/symits2


#pull the Defining Intragenomic Variants (DIVs) fasta file from the symportal website 

#wget https://symportal.org/static/resources/SymPortal_unique_DIVs.tar.gz --no-check-certificate

#tar -xvf SymPortal_unique_DIVs.tar.gz

#create a blast database 

#makeblastdb -in SymPortal_unique_DIVs.fasta  -dbtype nucl

#blast all reads.fq files against the ITS2 database 

DB="SymPortal_unique_DIVs.fasta"

# Function to convert FASTQ to FASTA
#convert_fastq_to_fasta() {
#    awk 'NR%4==1{sub(/^@/,">"); print} NR%4==2' "$1" > "${1%.fq}.fasta"
#}

#export -f convert_fastq_to_fasta

# Convert all .fq files to .fasta in parallel
#ls /data/putnamlab/tconn/apul_reseq/trimmed/*.fq | parallel convert_fastq_to_fasta {}


ls /data/putnamlab/tconn/apul_reseq/trimmed/*.fasta | parallel -j 20 blastn -task blastn -query {} -db "$DB" \
-num_threads 20 \
-outfmt 6 \
-evalue 1e-10 \
-out {.}.its2.blast.out
