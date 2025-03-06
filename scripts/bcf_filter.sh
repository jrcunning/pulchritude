#!/bin/bash
#SBATCH --job-name=bcf_filt
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --account=putnamlab
#SBATCH -o bcf_filt.out
#SBATCH -e bcf_filt.error


source ~/.bashrc
cd /data/putnamlab/tconn/apul_reseq/trimmed/bam

conda activate sambamba 

#filter sites in the vcf 

vcftools --vcf apulchra_moorea_2022.vcf --min-meanDP 5 \
--maf 0.5 \
--remove-indels \
--min-alleles 2 \
--max-alleles 2 \
--minQ 20 \
--max-missing 0.9 \
--out apulchra_moorea_2022_thin.vcf
