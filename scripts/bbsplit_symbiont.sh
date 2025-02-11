#!/bin/bash
#SBATCH --job-name=bbsplit
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=100gb
#SBATCH --time=24:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

source ~/.bashrc
conda  activate bbmap
module load parallel/20220722-GCCcore-11.3.0 

cd  /data/putnamlab/tconn/apul_reseq/trimmed 

for i in $(ls -1 *.trimmed.R1.fastq.gz | sed 's/\.trimmed.R1.fastq.gz//'); do

#skip the bam file if there are already symbiont files 
	sym_file=${i}_a.fq
	if [ -f "$sym_file" ]; then 
	echo "Skipping ${i}, symbiont file already exists"
	continue 
	fi

bbsplit.sh in1=$i\.trimmed.R1.fastq.gz in2=$i\.trimmed.R2.fastq.gz \
ref_a=/data/putnamlab/REFS/Sym_A1/Smic.genome.scaffold.final.fasta ref_b=/data/putnamlab/tconn/refs/Bminutum_genomic.fna ref_c=/data/putnamlab/REFS/Sym_C1/SymbC1.Genome.Scaffolds.fasta ref_d=/data/putnamlab/jillashey/Apul_Genome/dbs/102_symbd_genome_scaffold.fa basename=${i}_%.fq  --nodisk


done 
