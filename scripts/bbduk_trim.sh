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
cd  /data/putnamlab/tconn/apul_reseq

mkdir trimmed 

DIR=/data/putnamlab/tconn/apul_reseq/H2WLFDSXF/8740
OUT=/data/putnamlab/tconn/apul_reseq/trimmed

parallel --plus echo {%R._001.fastq.gz} ::: *.fastq.gz | parallel 'bbduk.sh \
in1=$DIR/.fq.gz \
in2=$DIR/Unknown_BK461-020033_2.fq.gz \
out1=$OUT/20033_1trimmed.fq.gz \
out2=$OUT/20033_2trimmed.fq.gz \
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
minlen=50 
