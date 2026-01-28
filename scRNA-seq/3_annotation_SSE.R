library(Seurat)
library(tibble)
library(dplyr)
library(patchwork)
library(celda)

SSE1_filtered <- readRDS("SSE1_4_seurat_umap_noAR_sse7.rds")

# Neftel subtypes
library(scalop)
load('Signatures_GBM.rda')
Signatures_GBM$G1S <- NULL
Signatures_GBM$G2M <- NULL

neo_score<-sigScores(m = as.matrix(SSE1_filtered[['RNA']]$counts),sigs = Signatures_GBM)

neo_score[1:4,1:4]
write.csv(neo_score,file = "Neftel_score_sse7.csv")
four_state<-as_four_state_gbm(as.matrix(neo_score))
four_state[1:4,1:4]

four_state$Neftel_subtype <- apply(four_state, 1, function(t) colnames(four_state)[which.max(t)])
head(four_state)

write.csv(four_state,file = "Neftel_subtype_sse7.csv")

# identify 'hybrids'
four_state <- four_state %>% 
  rownames_to_column(var = "ID") %>%
  rowwise() %>%
  mutate(
    # Find the scores of the top 3 modules for each cell
    top_modules = list(sort(c(OPC, NPC, AC, MES), decreasing = TRUE)[1:3]),
    
    # Find the names of the top 3 modules
    top_module_names = list(names(sort(c(OPC = OPC, NPC = NPC, AC = AC, MES = MES), decreasing = TRUE)[1:3]))
  ) %>%
  ungroup()

percentile_10 <- four_state %>%
  mutate(top_module = sapply(top_module_names, `[[`, 1)) %>%  # Extract top module for each cell
  group_by(top_module) %>%
  summarise(
    OPC = quantile(OPC[top_module == "OPC"], 0.1, na.rm = TRUE),
    NPC = quantile(NPC[top_module == "NPC"], 0.1, na.rm = TRUE),
    AC  = quantile(AC[top_module == "AC"], 0.1, na.rm = TRUE),
    MES = quantile(MES[top_module == "MES"], 0.1, na.rm = TRUE)
  )

four_state[i,9] <- NA
four_state[i,10] <- NA
for (i in 1:length(rownames(four_state))) {
  top_module_name = unlist(four_state[i,"top_module_names"])[1]
  second_module_score = unlist(four_state[i,"top_modules"])[2] 
  second_module_name  = unlist(four_state[i,"top_module_names"])[2] 
  third_module_score  = unlist(four_state[i,"top_modules"])[3]
  second_module_percentile = as.numeric(na.omit(percentile_10[[second_module_name]]))
  
  four_state[i,9] <- ifelse(second_module_score>1 &
                              (second_module_score-third_module_score)>0.3 &
                              second_module_score>second_module_percentile,'hybrids',four_state[i,"Neftel_subtype"])
  four_state[i,10] <- ifelse(four_state[i,9]=="hybrids",paste0(top_module_name,"_",second_module_name),NA)
}
colnames(four_state)[6] <- "Neftel_subtype_noHybrids"
colnames(four_state)[9] <- "Neftel_subtype"
colnames(four_state)[10] <- "hybrids"
four_state <- four_state[,-c(7,8)]
four_state <- as.data.frame(four_state)
four_state <- column_to_rownames(four_state,var = "ID")

table(four_state$Neftel_subtype)

write.csv(four_state,file = "Neftel_subtype_sse7_hybrids.csv")

SSE1_filtered <- AddMetaData(SSE1_filtered,four_state)



#steve pollard_imm subtype
library(AUCell)
library(GSEABase)
expr <- AUCell_buildRankings(SSE1_filtered[['RNA']]$counts,plotStats=F)

load("SPollard_imm.Rdata")
AUC_calculation <- AUCell_calcAUC(MESimm_SP,expr)
AUC_assignment <- AUCell_exploreThresholds(AUC_calculation,plotHist = T,assign=T)
AUC_SP <- AUC_calculation@assays@data@listData$AUC
AUC_SP <- as.data.frame(t(AUC_SP))

AUC_SP$SPollard_subtype <- apply(AUC_SP, 1, function(t) colnames(AUC_SP)[which.max(t)])
head(AUC_SP)

SSE1_filtered <- AddMetaData(SSE1_filtered,AUC_SP)
save(AUC_calculation,AUC_assignment,AUC_SP,file = "AUC_SP.RData")



#development injury
load("Developmental_injury.Rdata")
AUC_calculation <- AUCell_calcAUC(Developmental_injury,expr)
AUC_assignment <- AUCell_exploreThresholds(AUC_calculation,plotHist = T,assign=T)
AUC_DI <- AUC_calculation@assays@data@listData$AUC
AUC_DI <- as.data.frame(t(AUC_DI))

AUC_DI$Dirks_subtype <- apply(AUC_DI, 1, function(t) colnames(AUC_DI)[which.max(t)])
head(AUC_DI)

SSE1_filtered <- AddMetaData(SSE1_filtered,AUC_DI)
save(AUC_calculation,AUC_assignment,AUC_DI,file = "AUC_DI.RData")



#cell cycle phase
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
SSE1_filtered=CellCycleScoring(object = SSE1_filtered, 
                               s.features = s.genes, 
                               g2m.features = g2m.genes, 
                               set.ident = TRUE)



#cell cycle ccAFv2
library(keras)
library(ccAFv2)
library(reticulate)
use_condaenv("/exports/cmvm/eddie/scs/groups/spollar2-PollardLab/Zeyu/anaconda/envs/ccAFv2_R")

SSE1_filtered <- PredictCellCycle(SSE1_filtered,gene_id='symbol')
ccAF_phase <- SSE1_filtered@meta.data[,c("G1","G2.M","Late.G1","M.Early.G1","Neural.G0","S","S.G2","ccAFv2")]
  
write.csv(ccAF_phase,file = "ccAF_results_combine_sse7.csv")

SSE1_filtered <- AddMetaData(SSE1_filtered,ccAF_phase)

write.csv(SSE1_filtered@meta.data,file = "meta_annotation_sse7.csv")

saveRDS(SSE1_filtered,file = "SSE1_5_seurat_final_sse7.rds")
