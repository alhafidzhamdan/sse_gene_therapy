#!/bin/bash

#$ -N consensus
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=2G
#$ -l h_rt=00:05:00

### First generate bed files containing overlapping SOX2 and SOX9 peaks for each of the 7 samples:
### and sort by chromosome and start position:

for file in ChIP-seq/processed_data/*_Sox2_*narrowPeak; do
    sample=$(basename $file | cut -d'_' -f1)
    echo $sample
    bedtools intersect -a $file -b ChIP-seq/processed_data/${sample}_Sox9_peaks.narrowPeak | sort -k1,1 -k2,2n > ChIP-seq/processed_data/${sample}_Sox2_Sox9_overlap.bed
done

### Using bedtools multiintersect to find consensus overlapping SOX2 and SOX9 peaks across replicates:
bedtools multiinter -header -i ChIP-seq/processed_data/*_Sox2_Sox9_overlap.bed  > ChIP-seq/processed_data/Merged_rep_overlap.bed


