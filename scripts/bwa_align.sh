#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

source ~/.bashrc

##load bwa
conda  activate sambamba


cd ~/scratch

bwa mem -t 20  ~/scratch/ofav_comparisons/um_ofav_v1_softmask.fa   ~/scratch/{} ~/scratch/{} |samtools view -S -b  -o - - | samtools sort -o {}.bam



