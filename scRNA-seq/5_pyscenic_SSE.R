library(Seurat)
library(tibble)
library(dplyr)
library(SCopeLoomR)

SSE1_filtered <- readRDS("SSE1_5_seurat_final_sse7.rds")
SSE1_filtered <- subset(SSE1_filtered,Condition_mCherry_strict %in% c("Control_No","SSE7_No","SSE7_Yes"))

# scenic by library
for (i in unique(SSE1_filtered$library)) {
  a <- subset(SSE1_filtered,library == i)
  #build_loom
  path <- paste0("./scenic/",i)
  if (!dir.exists(path)) {
    dir.create(path)
  }
  
  build_loom(file.name=paste0(path,"/data.loom"),dgem=a[['RNA']]$counts)
}

# scenic by cell line
for (i in c("E17","E20","E21","E28","E31","E34","E43","E55")) {
  a <- readRDS(paste0("SSE1_6_seurat_",i,"_sse7.rds"))
  #build_loom
  path <- paste0("./scenic/",i)
  if (!dir.exists(path)) {
    dir.create(path)
  }
  
  build_loom(file.name=paste0(path,"/data.loom"),dgem=a[['RNA']]$counts)
}


## run this part with array job 
## command on eddie
qsub -t 1-15 SSE_array_lib.sh

qsub -t 1-8 SSE_array_line.sh



