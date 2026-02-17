#!/bin/bash

## Example bash script to submit a job to compute the matrix and plot a heatmap of signal over peaks using deepTools. To run this script, do: qsub -t 1-n submit_deeptools_heatmap.sh CONFIG IDS
## CONFIG is the path to the file scripts/config.sh which contains environment variables set to commonly used paths and files in the script
## IDS is a list of sample ids, one per line.

#$ -N plotHeatmapProfile
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=12G
#$ -l h_rt=1:00:00
#$ -pe sharedmem 12
    
CONFIG=$1
IDS=$2

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1`

### Needs deepTools which needs python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

### Plot SOX9 signal in SOX2 peaks (and vice versa) for each sample. Here we plot SOX9 signal in SOX2 peaks. To plot SOX2 signal in SOX9 peaks, just swap the region bed and the signal track variables below.
BASE_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF
OUTPUT_NAME=${SAMPLE_ID}_Sox9_signal_over_Sox2_peaks
OUTPUT_DIR=${BASE_DIR}/qc/deeptools/computeMatrix/$OUTPUT_NAME
REGION_BED=${BASE_DIR}/downstream/macs2/merged/${SAMPLE_ID}_Sox2_peaks.narrowPeak
SIGNAL_TRACK=${BASE_DIR}/downstream/coverage/bigwig/${SAMPLE_ID}_Sox9.normalised.bw

### Run:
if [[ ! -d $OUTPUT_DIR ]]; then

    mkdir -p $OUTPUT_DIR

fi

computeMatrix reference-point \
    --referencePoint center \
    -R $REGION_BED \
    -S $SIGNAL_TRACK \
    -b 5000 -a 5000 \
    --binSize 10 \
    --skipZeros --missingDataAsZero \
    --blackListFileName $BLACKLISTED_PEAKS \
    -p 32 \
    --scale 10 \
    --smartLabels \
    -o ${OUTPUT_DIR}/${OUTPUT_NAME}.matrix.gz \
    --outFileNameMatrix ${OUTPUT_DIR}/${OUTPUT_NAME}.matrix.tab \
    --outFileSortedRegions ${OUTPUT_DIR}/${OUTPUT_NAME}.matrix.bed 

if [[ -f ${OUTPUT_DIR}/${OUTPUT_NAME}.matrix.gz ]]; then

    plotHeatmap \
        -m ${OUTPUT_DIR}/${OUTPUT_NAME}.matrix.gz \
        -out ${OUTPUT_DIR}/${OUTPUT_NAME}.heatmap.pdf \
        --plotFileFormat pdf \
        --dpi 300 \
        --zMin -3 --zMax 3 \
        --legendLocation none \
        --verbose 

fi

    
    
    
    
    
    
