#!/bin/bash
#SBATCH --job-name=ngsrelate
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=200gb
#SBATCH --time=144:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

source ~/.bashrc
PATH=$PATH
REF=/data/putnamlab/tconn/annotate_results/final_files/Acropora_pulchra_v1.1.fa.masked 
conda activate angsd 
module load NgsRelate 


cd /data/putnamlab/tconn/apul_reseq/trimmed/bam

##extract frequency column from the previously estimated allele frequency file 

zcat genopulchra.mafs.gz | cut -f5 | sed 1d >pulchrafreq


##run NgsRelate to estimate relatedness coefficients between samples of Acropora puclrha 

ngsRelate -g genopulchra.glf.gz -n 206 -f pulchrafreq -O pulchra_ngsrelate  -p 20
