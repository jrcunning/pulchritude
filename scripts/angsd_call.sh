#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100gb
#SBATCH --time=96:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error


#variant calling with ANGSD 

source ~/.bashrc
PATH=$PATH
REF=/data/putnamlab/tconn/annotate_results/final_files/Acropora_pulchra_v1.1.fa.masked 
conda activate angsd 

#load angsd 
cd /data/putnamlab/tconn/apul_reseq/trimmed/bam

ls *_nodup.bam_final.bam > bamlist

angsd -GL 1 \ #use samtools calling method -- GATK not as good for non-model organisms
 -out genopulchra \
 -nthreads 10 \
 -doGLF 3 \ #output beagle file (what do we need for PCAangsd/SNP calling)
 -bam bamlist \
 -domajorminor 1 \
 -snp_pval 1e-6 \ #p-value for calling a snp 
 -domaf 1 \
 -minmaf 0.05 # minimum minor allele frequency  

