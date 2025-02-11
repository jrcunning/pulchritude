#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=100gb
#SBATCH --time=96:00:00
#SBATCH --account=putnamlab
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.error

source ~/.bashrc
PATH=$PATH
##load bwa
conda  activate sambamba
module load picard 

#bwa index /data/putnamlab/tconn/annotate_results/final_files/Acropora_pulchra_v1.1.fa.masked 

cd /data/putnamlab/tconn/apul_reseq/trimmed/bam

#loop through files 
for i in *_L003_R1_001.fastq.gz.bam; do

#run picard to change read groups 
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$i O=${i}_rg.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=mo22 RGSM=$i 

#remove duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${i}_rg.bam O=${i}_nodup.bam M=$i_marked_dup.txt REMOVE_SEQUENCING_DUPLICATES=true

done

