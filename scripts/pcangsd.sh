#!/bin/bash
#SBATCH --job-name=pcangsd
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
module load PCAngsd
#load angsd 
cd /data/putnamlab/tconn/apul_reseq/trimmed/bam

#angsd -GL 1 -out genopulchra -nthreads 10 -doGLF 3 -bam bamlist -domajorminor 1 -snp_pval 1e-6 -doMaf 1 -minmaf 0.05  
#filter sites before running pcangsd, keeping 

#run pcangsd, have the program select the number of eigenvectors, ignore selection flag for now while it is a bunch of similarly related individuals

pcangsd --beagle genopulchra.beagle.gz --out pulchraout --threads 20 --iter 500


