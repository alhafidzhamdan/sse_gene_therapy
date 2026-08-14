# Synthetic super-enhancers enable precision viral immunotherapy
This repository contains steps taken to analyse single cell RNA-seq and SOX2 and SOX9 ChIP-seq data in GBM stem-like cells (GSCs), as used in Koeber and Matjusaitis et al, Nature 2026.

## scRNA-seq

Steps:
8 human GSC lines (E17, E20, E21, E28, E31, E34, E43 and E55) were transduced with AAV-SSE7-mCherry. 3 days post-transduction, 10% of the cells were collected for flow cytometry to assess transduction efficiency, and the remaining 90% were processed for single-cell RNA-sequencing library preparation using the PARSE Evercode v2 Cell Fixation Kit and the PARSE Evercode WT Mega v2 Kit. Paired-end reads were aligned to the hg38 genome using the standard PARSE Bioscience pipeline, split-pipe (version 1.1.2), with default parameters. Three custom genes: mCherry, HSV-TK, and bGH polyA; were added to the reference genome to enable detection of activated SSE-7 cells.


Data: 
- Raw fastqs files can be obtained from the ENA repository under the study number PRJEB81816


## ChIP-seq

Steps:
Paired-end reads were aligned to the hg38 genome using BWA, filtering out poor quality (MAPQ <10), duplicates, mitochondrial genome and blacklisted regions (ENCSR636HFF). Technical replicates were merged, and peaks were called using MACS2 with default settings. Consensus SOX2 and SOX9 peak sets were derived by taking the overlap of peaks occurring in 5 out of 7 GSCs for each TF set. The overlap significance between the consensus SOX2 and SOX9 peak sets was determined using a circular permutation test in regioneR (ntimes = 10000). To perform overlapping analysis with GSC SEs, we downloaded raw fastq files from publicly available H3K27ac datasets on GSCs (GSE119834, GSE74529, GSE121601, and GSE92458) and called SEs using ROSE with default settings. A consensus set of GSC SEs were derived by considering SEs that occurred in at least two GSCs. The overlap significance between the consensus co-bound SOX2/SOX9 sites and consensus GSC SEs was performed using a circular permutation test as above. GO term association analysis was performed using GREAT. To find centrally enriched de novo motifs at SOX2 and SOX9 peaks within co-bound SOX2/SOX9 enhancers, and to identify the spacing between the most significant de novo motifs from each peak set, we used CentriMo and SpaMo, from the MEME-ChIP suite of tools.

Data: 
- Raw (fastqs) and processed data (bigwigs and merged peak files) can be obtained from the ENA repository under the study number PRJEB107008
- Merged peak files (and consensus peak sets) can be also obtained within this repository under ChIP-seq/processed_data
- SOX2 and SOX9 ChIP-seq coverage tracks are hosted at https://genome.ucsc.edu/s/alhafidzhamdan/co_bound_SOX2_SOX9_peaks


## Citation

    URL: https://www.nature.com/articles/s41586-026-10329-6
    DOI: https://doi.org/10.1038/s41586-026-10329-6


