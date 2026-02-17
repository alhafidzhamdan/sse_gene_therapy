#!/bin/bash

## Example bash script to call peaks using MACS2 on merged replicates. To run this script, do: qsub -t 1-n submit_macs2.sh CONFIG IDS BAM_DIR
## CONFIG is the path to the file scripts/config.sh which contains environment variables set to commonly used paths and files in the script.
## IDS is a list of sample ids, one per line.
## BAM_DIR is the path to the directory containing the final bams for each sample and their corresponding input bams.

#$ -N macs2
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
#$ -pe sharedmem 8

CONFIG=$1
IDS=$2
BAM_DIR=$3

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1`

### Merge replicate ChIP bams:
samtools merge -@ 10 ${BAM_DIR}/${SAMPLE_ID}_merged.final.bam ${BAM_DIR}/${SAMPLE_ID}_Rep1.final.bam ${BAM_DIR}/${SAMPLE_ID}_Rep2.final.bam
samtools index ${BAM_DIR}/${SAMPLE_ID}_merged.final.bam

### Call peaks using MACS2 with the merged ChIP bam and the singular input bam as control.
macs2 callpeak \
    -t ${BAM_DIR}/${SAMPLE_ID}_merged.final.bam \
    -c ${BAM_DIR}/${SAMPLE_ID}_Input_Rep1.final.bam \
    -n ${SAMPLE_ID} \
    -B \
    -p 1e-5








