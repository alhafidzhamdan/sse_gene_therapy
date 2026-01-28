library(Seurat)
library(tibble)
library(dplyr)

SSE1_filtered <- readRDS("SSE1_5_seurat_final_sse7.rds")
SSE1_filtered <- subset(SSE1_filtered,Condition_mCherry_strict %in% c("Control_No","SSE7_No","SSE7_Yes"))

for (i in unique(SSE1_filtered$cell_line)) {
  a <- subset(SSE1_filtered,cell_line == i)
 
  a <- NormalizeData(a, normalization.method = "LogNormalize")
  a <- FindVariableFeatures(a, selection.method = "vst", 
                            nfeatures = 2000)
  a <- ScaleData(a, vars.to.regress = c("percent.mito"))
  a <- RunPCA(a, seed.use = 240508)
  
  a <- FindNeighbors(a, dims = 1:15) 
  a <- FindClusters(object = a, verbose = T, resolution = 1)
  a <- RunUMAP(object = a, dims = 1:15,
               n.neighbors = 50,min.dist = .25,
               n.epochs = 300,seed.use = 240508)
  
  saveRDS(a,file = paste0("SSE1_6_seurat_",i,"_sse7.rds"))
}

