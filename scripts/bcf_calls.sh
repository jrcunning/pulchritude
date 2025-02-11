#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100gb
#SBATCH --time=96:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

##script to call variants using bcftools

source ~/.bashrc
PATH=$PATH
REF=/data/putnamlab/tconn/annotate_results/final_files/Acropora_pulchra_v1.1.fa.masked 

#activate bcftools conda environment 
conda  activate sambamba
#change directory 
cd /data/putnamlab/tconn/apul_reseq/trimmed/bam
#create list of bam files to call from 

ls *_nodup.bam_final.bam > bamlist

#pipe mpileup to calls 

bcftools mpileup -f $REF -b bamlist | bcftools call -mv -o apulchra_moorea_2022.vcf


