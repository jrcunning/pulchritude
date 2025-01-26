#!/bin/bash
#SBATCH --job-name=bbsplit
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

source ~/.bashrc
conda  activate bbmap
module load parallel/20220722-GCCcore-11.3.0 

cd  /data/putnamlab/tconn/apul_reseq/trimmed 

parallel -j 20 bbsplit.sh \
in1={} \
in2={=s/R1/R2/=}\
ref_a=/data/putnamlab/REFS/Sym_A1/Smic.genome.scaffold.final.fasta ref_b=/data/putnamlab/tconn/refs/Bminutum_genomic.fna ref_c=/data/putnamlab/REFS/Sym_C1/SymbC1.Genome.Scaffolds.fasta ref_d=/data/putnamlab/jillashey/Apul_Genome/dbs/102_symbd_genome_scaffold.fa out_a={}_a.fq out_b={}_b.fq out_c={}_c.fq out_d={}_d.fq ::: *.trimmed.R1.fastq.gz
