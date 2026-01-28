library(Seurat)
library(tibble)
library(dplyr)
library(DoubletFinder)
library(celda)

SSE1_filtered <- readRDS("SSE1_1_seurat_filtered_sse7.rds")

# UMAP: standard or SCT
# standard
SSE1_filtered <- NormalizeData(SSE1_filtered, normalization.method = "LogNormalize")
SSE1_filtered <- FindVariableFeatures(SSE1_filtered, selection.method = "vst", 
                                      nfeatures = 2000)
SSE1_filtered <- ScaleData(SSE1_filtered, vars.to.regress = c("percent.mito"))

# SCT
#SSE1_filtered <- SCTransform(SSE1_filtered, method = "glmGamPoi",vars.to.regress = "percent.mito",seed.use = 240508)

# after standard or SCT
SSE1_filtered <- RunPCA(SSE1_filtered,seed.use = 240508)
p1 <- ElbowPlot(SSE1_filtered, ndims=50, reduction="pca") 
ggsave("QC_elbow_sketch.pdf",p1)

SSE1_filtered <- FindNeighbors(SSE1_filtered, dims = 1:15) 
SSE1_filtered <- FindClusters(object = SSE1_filtered, verbose = T, resolution = 1)
SSE1_filtered <- RunUMAP(object = SSE1_filtered, dims = 1:15,
                         n.neighbors = 50,min.dist = .25,
                         n.epochs = 300,seed.use = 240508)

p1 <- DimPlot(SSE1_filtered,
              reduction = "umap", 
              group.by = 'seurat_clusters',
              pt.size = 0.5,label = T)
ggsave("umap_seurat_clusters_sse7.pdf",p1)

saveRDS(SSE1_filtered,file = "SSE1_2_seurat_umap_sse7.rds")

# doublets
doublets_final <- data.frame()
for (i in unique(SSE1_filtered$library)) {
  a <- subset(SSE1_filtered,library==i)
  
  a <- NormalizeData(a, normalization.method = "LogNormalize")
  a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)
  a <- ScaleData(a)
  a <- RunPCA(a,seed.use = 240508)
  a <- RunUMAP(object = a, dims = 1:15,
               n.neighbors = 50,min.dist = .25,
               n.epochs = 300,seed.use = 240508)
  
  # doublets
  ## best pk
  sweep.res.list <- paramSweep(a, PCs = 1:15, sct = FALSE)
  
  #non ground truth
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  DoubletRate = 0.03    
  homotypic.prop <- modelHomotypic(a$sample)  
  nExp_poi <- round(DoubletRate*ncol(a)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  
  a <- doubletFinder(a, PCs = 1:15, pN = 0.25, pK = pK_bcmvn, 
                     nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  new_pANN <- paste0("pANN_0.25_",pK_bcmvn,"_",nExp_poi)
  a <- doubletFinder(a, PCs = 1:15, pN = 0.25, pK = pK_bcmvn, 
                     nExp = nExp_poi.adj, reuse.pANN = new_pANN, 
                     sct = FALSE)

  
  b <- a@meta.data[,30:33]
  colnames(b)[2:4] <- c("parameters","class_heter","class_homo")
  if (length(doublets_final)==0) {
    doublets_final <- b
  } else {doublets_final <- rbind(doublets_final,b)}
}

write.csv(doublets_final,file = "meta_doubletsfinder_nonGT.csv")

meta_doublets <- read.csv("./meta_info/meta_doubletsfinder_nonGT_sse7.csv",row.names = 1)
SSE1_filtered <- AddMetaData(SSE1_filtered,meta_doublets)

table(SSE1_filtered$class_homo,SSE1_filtered$seurat_clusters)

SSE1_filtered <- subset(SSE1_filtered,class_homo == "Singlet")
clusters_keep <- setdiff(unique(SSE1_filtered$seurat_clusters), c(21,23:26,28,29,32)) # remove cluster with over 20% doublets
SSE1_filtered <- subset(SSE1_filtered,seurat_clusters %in% clusters_keep) 

SSE1_singlet <- colnames(SSE1_filtered)
write.table(SSE1_singlet,file = "SSE1_Singlet_cellname_sse7.txt")

saveRDS(SSE1_filtered,file = "SSE1_3_seurat_umap_noD_sse7.rds")

decontX_results <- decontX(SSE1_filtered[['RNA']]$counts)
SSE1_filtered$Contamination =decontX_results$contamination

SSE1_filtered = SSE1_filtered[,SSE1_filtered$Contamination < 0.2]

saveRDS(SSE1_filtered,file = "SSE1_4_seurat_umap_noAR_sse7.rds")
