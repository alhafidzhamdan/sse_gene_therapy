# Synthetic super-enhancers enable selective expression of anti-cancer payloads for viral gene therapy
This repository contains steps taken to analyse single cell RNA-seq data and SOX2 and SOX9 ChIP-seq data on GBM stem-like cells (GSCs), used in Koeber and Matjusaitis et al, Nature 2026.

## scRNA-seq



## ChIP-seq

Paired-end reads were aligned to the hg38 genome using BWA, filtering out poor quality (MAPQ <10), duplicates, mitochondrial genome and blacklisted regions (ENCSR636HFF). Technical replicates were merged, and peaks were called using MACS2 with default settings. Consensus SOX2 and SOX9 peak sets were derived by taking the overlap of peaks occurring in 5 out of 7 GSCs for each TF set. The overlap significance between the consensus SOX2 and SOX9 peak sets was determined using a circular permutation test in regioneR (ntimes = 10000). Genomic regions near SOX2 and SOX9 peaks were annotated using HOMER. 

To perform overlapping analysis with GSC SEs, we downloaded raw fastq files from publicly available H3K27ac datasets on GSCs (GSE119834, GSE74529, GSE121601, and GSE92458) and called SEs using ROSE with default settings. A consensus set of GSC SEs were derived by considering SEs that occurred in at least two GSCs. The overlap significance between the consensus co-bound SOX2/SOX9 sites and consensus GSC SEs was performed using a circular permutation test as above. GO term association analysis was performed using GREAT. To find centrally enriched de novo motifs at SOX2 and SOX9 peaks within co-bound SOX2/SOX9 enhancers, and to identify the spacing between the most significant de novo motifs from each peak set, we used CentriMo and SpaMo, from the MEME-ChIP suite of tools (Ma et al., 2014). 



Citation

    URL: Pending
    DOI: Pending


