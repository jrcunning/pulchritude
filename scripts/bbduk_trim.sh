#!/bin/bash
#SBATCH --job-name=bbduktrim
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

cd  /data/putnamlab/tconn/apul_reseq/H2WLFDSXF/8740

parallel -j 20 bbduk.sh \
in1={} \
in2={=s/R1/R2/=} \
out1=/data/putnamlab/tconn/apul_reseq/trimmed/{}.trimmed.R1.fastq.gz \
out2=/data/putnamlab/tconn/apul_reseq/trimmed/{}.trimmed.R2.fastq.gz \
ref=/home/trinity.conn/.conda/envs/bbmap/opt/bbmap-39.13-1/resources/adapters.fa \
ktrim=r \
k=23 \
mink=11 \
hdist=1 \
tpe \
tbo \
qtrim=r \
tossbrokenreads=t \
trimq=20 \
maq=20 \
minlen=50 ::: *_L003_R1_001.fastq.gz
