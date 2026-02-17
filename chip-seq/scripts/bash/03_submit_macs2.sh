#!/bin/bash

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








