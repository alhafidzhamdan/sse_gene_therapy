library(Seurat)
library(BPCells)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scater)

DGE_folder <- "./Pollard_lab_SSE/"

mat <- ReadParseBio(DGE_folder)

# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(DGE_folder, "cell_metadata.csv"))

meta <- readxl::read_xlsx(paste0(DGE_folder, "meta.xlsx"))
rownames(meta) <- paste("SPZW",meta$SampleID,sep = "")
meta <- rownames_to_column(meta,var = "sample") %>% .[,-2]
colnames(meta)[2] <- c("cell_line")
colnames(meta)[5] <- c("mCherry_flow")

cell_meta <- merge(cell_meta,meta,by="sample") %>% column_to_rownames(var = "bc_wells")
cell_meta$library <- paste0("library",as.numeric(substr(rownames(cell_meta),12,13)))

# Create object
SSE1 <- CreateSeuratObject(mat, meta.data = cell_meta)
SSE1_SSE7 <- subset(SSE1,Conditions %in% c("Control","SSE7"))
saveRDS(SSE1_SSE7,"SSE1_0_seurat_sse7.rds")


## data QC
C<-GetAssayData(object = SSE1_SSE7, slot = "counts")

# mito
mito_genes=rownames(SSE1_SSE7)[grep("^MT-", rownames(SSE1_SSE7))]
percent.mito <- Matrix::colSums(C[mito_genes,])/Matrix::colSums(C)*100
SSE1_SSE7 <- AddMetaData(SSE1_SSE7,percent.mito,col.name = "percent.mito")

# ribo
ribo_genes=rownames(SSE1_SSE7)[grep("^RP[sl]", rownames(SSE1_SSE7),ignore.case = T)]
percent.ribo <- Matrix::colSums(C[ribo_genes,])/Matrix::colSums(C)*100
SSE1_SSE7 <- AddMetaData(SSE1_SSE7,percent.ribo,col.name = "percent.ribo")

# hb
hb_genes=rownames(SSE1_SSE7)[grep("^HB[^(p)]", rownames(SSE1_SSE7))]
percent.hb <- Matrix::colSums(C[hb_genes,])/Matrix::colSums(C)*100
SSE1_SSE7 <- AddMetaData(SSE1_SSE7,percent.hb,col.name = "percent.hb")


summary(SSE1_SSE7$nCount_RNA)
summary(SSE1_SSE7$nFeature_RNA)
summary(SSE1_SSE7$log10GenesPerUMI)
summary(SSE1_SSE7$percent.mito)
summary(SSE1_SSE7$percent.ribo)
summary(SSE1_SSE7$percent.hb)


p1 <- FeatureScatter(SSE1_SSE7, raster = T,
                     feature1 = "nCount_RNA", 
                     feature2 = "nFeature_RNA",
                     group.by = "library")
p2 <- FeatureScatter(SSE1_SSE7, raster = T, 
                     feature1 = "nCount_RNA", 
                     feature2 = "percent.mito",
                     group.by = "library")
p3 <- FeatureScatter(SSE1_SSE7, raster = T, 
                     feature1 = "nCount_RNA", 
                     feature2 = "percent.ribo",
                     group.by = "library")
p4 <- FeatureScatter(SSE1_SSE7, raster = T, 
                     feature1 = "nCount_RNA", 
                     feature2 = "percent.hb",
                     group.by = "library")

ggplot(SSE1@meta.data) +
  ggridges::geom_density_ridges(aes(x = percent.mito, y = cell_line, 
                                    fill = Conditions),alpha=0.5)

SSE1_filtered_SSE7 <- subset(SSE1_SSE7,percent.mito < 20 & 
                          percent.hb < 10 &
                          percent.ribo < 30 &
                          nFeature_RNA > 300 &
                          #nFeature_RNA < 7500 & #depends on cor between nFeature and nCount
                          nCount_RNA > 500 & # UMI
                          nCount_RNA < 20000)


counts <- GetAssayData(object = SSE1_filtered_SSE7, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) > 5

SSE1_filtered_SSE7 <- SSE1_filtered_SSE7[keep_genes,]


SSE1_filtered_SSE7$group <- paste(SSE1_filtered_SSE7$cell_line,SSE1_filtered_SSE7$Conditions,sep = "_")
SSE1_filtered_SSE7$line_passage <- paste(SSE1_filtered_SSE7$cell_line,SSE1_filtered_SSE7$Passage,sep = "_")
SSE1_filtered_SSE7$group_passage <- paste(SSE1_filtered_SSE7$group,SSE1_filtered_SSE7$Passage,sep = "_")

SSE1_filtered_SSE7$mCherry_count <- ifelse(SSE1_filtered_SSE7[['RNA']]$counts["mCherry",]>0 | 
                                        SSE1_filtered_SSE7[['RNA']]$counts["HSV-TK",]>0 | 
                                        SSE1_filtered_SSE7[['RNA']]$counts["bGHployA",]>0, "Yes","No")
SSE1_filtered_SSE7$Condition_mCherry <- paste(SSE1_filtered_SSE7$Conditions,SSE1_filtered_SSE7$mCherry_count,sep = "_")

SSE1_filtered_SSE7$mCherry_count_strict <- ifelse(SSE1_filtered_SSE7[['RNA']]$counts["mCherry",]>1 | 
                                               SSE1_filtered_SSE7[['RNA']]$counts["HSV-TK",]>1 | 
                                               SSE1_filtered_SSE7[['RNA']]$counts["bGHployA",]>0, "Yes","No")

SSE1_filtered_SSE7$Condition_mCherry_strict <- paste(SSE1_filtered_SSE7$Conditions,SSE1_filtered_SSE7$mCherry_count_strict,sep = "_")

SSE1_filtered_SSE7$Condition_mCherry_strict <- ifelse(SSE1_filtered_SSE7$mCherry_count_strict =="No" &
                                                 SSE1_filtered_SSE7$Conditions != "Control" &
                                                 SSE1_filtered_SSE7$mCherry_count =="Yes" , 
                                                 "False_Neg",SSE1_filtered_SSE7$Condition_mCherry_strict)

SSE1_filtered_SSE7$Condition_mCherry_strict[which(SSE1_filtered_SSE7$mCherry_count=="Yes" &
                                                    SSE1_filtered_SSE7$Conditions == "Control")] <- "False_Pos"


saveRDS(SSE1_filtered_SSE7,"SSE1_1_seurat_filtered_SSE7.rds")

