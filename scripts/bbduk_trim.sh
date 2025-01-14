#!/bin/bash
#SBATCH --job-name=bbduktrim
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --output=trim.out
#SBATCH --output=trim.error

source ~/.bashrc
conda  activate bbmap
module load parallel/20220722-GCCcore-11.3.0 
cd  /data/putnamlab/tconn/apul_reseq/H2WLFDSXF/8740

DIR=/data/putnamlab/tconn/apul_reseq/H2WLFDSXF/8740
OUT=/data/putnamlab/tconn/apul_reseq/trimmed

parallel -j 20 bbduk.sh \
in1={} \
in2={=s/R1/R2/=} \
out1=$OUT/{%.R1.fastq.gz}.trimmed.R1.fastq.gz \
out2=$OUT/{%.R1.fastq.gz}.trimmed.R2.fastq.gz \
ref=/data/putnamlab/tconn/apul_reseq/adapters.fa \
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
