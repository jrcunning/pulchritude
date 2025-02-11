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
PATH=$PATH
##load bwa
conda  activate sambamba

#bwa index /data/putnamlab/tconn/annotate_results/final_files/Acropora_pulchra_v1.1.fa.masked 

cd /data/putnamlab/tconn/apul_reseq/trimmed/bam

#loop through files 
for i in *_nodup.bam; do 

samtools view -q 20 -f 0x2 -b $i > ${i}_final.bam
samtools index ${i}_final.bam
samtools flagstat ${i}_final.bam > ${i}_stat.txt

done
