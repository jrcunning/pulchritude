#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100gb
#SBATCH --time=96:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

source ~/.bashrc

##load bwa
conda  activate sambamba

#bwa index /data/putnamlab/tconn/annotate_results/final_files/Acropora_pulchra_v1.1.fa.masked 

cd /data/putnamlab/tconn/apul_reseq/trimmed

#loop through files 
for i in $(ls -1 *.trimmed.R1.fastq.gz | sed 's/\.trimmed.R1.fastq.gz//'); do 
#skip  the files if a bam already exists 
	if [-f "${i}.bam" ]; then 
	echo "Skipping ${i}, BAM file already exists."
	continue 
	fi

#run bwa 
echo "Processing ${i}..."
bwa mem -t 20 /data/putnamlab/tconn/annotate_results/final_files/Acropora_pulchra_v1.1.fa.masked  $i\.trimmed.R1.fastq.gz  $i\.trimmed.R2.fastq.gz  |samtools view -S -b  -o - - | samtools sort -o $i.bam

done

