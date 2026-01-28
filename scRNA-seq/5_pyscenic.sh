#!/bin/sh

# set up module environment
. /etc/profile.d/modules.sh


module load anaconda
source activate pyscenic

cd ./scenic/library$1/
rm -r dask-worker-space

## expression matrix
f_ex_matrix_csv=data.loom

## referrence databases
f_db_names=../hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
f_motif_path=../motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
f_tfs=../allTFs_hg38.txt


pyscenic grn \
--output adj.tsv \
--method grnboost2 \
--seed 240607 \
--num_workers 32 \
${f_ex_matrix_csv} \
${f_tfs} 

pyscenic ctx adj.tsv \
${f_db_names} \
--annotations_fname ${f_motif_path} \
--expression_mtx_fname ${f_ex_matrix_csv} \
--output sce_regulon.gmt \
--mask_dropouts 

pyscenic aucell \
${f_ex_matrix_csv} \
sce_regulon.gmt \
--output sce_regulon_AUC.loom \
--seed 240607 
