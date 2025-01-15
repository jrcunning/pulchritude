#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

source ~/.bashrc

##load fastqc
conda activate fastqc 

cd /data/putnamlab/tconn/apul_reseq/trimmed

for i in *_L003_R1_001.fastq.gz.trimmed.*.fastq.gz; do 

fastqc $i 

done 


done

