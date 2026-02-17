#!/bin/bash

## Example bash script to align ChIP-seq reads to the reference genome using BWA. To run this script, do: qsub -t 1-n submit_bwa.sh CONFIG IDS READ_DIR BAM_DIR
## CONFIG is the path to the file scripts/config.sh which contains environment variables set to commonly used paths and files in the script.
## IDS is a list of sample ids, one per line.
## READ_DIR is the path to the directory containing the raw fastq files for each sample.
## BAM_DIR is the path to the directory where the aligned bams will be saved.

#$ -N bwa
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=90:00:00
#$ -pe sharedmem 16

CONFIG=$1
IDS=$2
READ_DIR=$3
BAM_DIR=$4

source $CONFIG

REPLICATE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 8`

# First run `bwa index $BWA_REF`
## Here we used hg38

# Align paired-end fastq reads 

cd $READ_DIR
echo "Aligning fastq file for ${REPLICATE_ID} ... "
bwa mem -M -t 20 $BWA_REF ${REPLICATE_ID}_R1.fastq.gz ${REPLICATE_ID}_R2.fastq.gz | samtools sort -@ 10 - -T ${REPLICATE_ID} -o $BAM_DIR/${REPLICATE_ID}.sorted.bam
samtools index $BAM_DIR/${REPLICATE_ID}.sorted.bam
samtools flagstat $BAM_DIR/${REPLICATE_ID}.sorted.bam > $BAM_DIR/${REPLICATE_ID}.stats.out
