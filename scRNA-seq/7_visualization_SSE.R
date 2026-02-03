library(Seurat)
library(BPCells)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scater)
library(RColorBrewer)
library(emmeans)
library(forcats)
source("custom_function.R")

# convert to bpcells
SSE1_filtered <- readRDS("SSE1_5_seurat_final_sse7.rds")

write_matrix_dir(mat = SSE1_filtered[['RNA']]$counts, dir = './bpcells_sse7')
counts.mat <- open_matrix_dir(dir = "./bpcells_sse7")

meta_annotation <- SSE1_filtered@meta.data

SSE1_filtered <- CreateSeuratObject(counts.mat, meta.data = meta_annotation)
saveRDS(SSE1_filtered,file = "SSE1_6_bpcells_final_sse7.rds")

# start analysis
SSE1_filtered <- readRDS("SSE1_6_bpcells_final_sse7.rds")
SSE1_filtered <- subset(SSE1_filtered,Condition_mCherry_strict %in% c("Control_No","SSE7_No","SSE7_Yes"))


# mCherry positive cells composition
library(scales)

# heterogeneity in ctr group
subtype_plot(SSE1_filtered,"ccAFv2","all", condition="Control_No")

for (i in c("Neftel_subtype","Phase","Dirks_subtype","SPollard_subtype","ccAFv2")) {
  p1 <- subtype_plot(SSE1_filtered,i,"all", condition="Control_No")
  ggsave(filename = paste0("Barplot_control_",i,".pdf"),p1)
}

# heterogeneity and SSE7 activation
cherry_proportion_plot(SSE1_filtered,"ccAFv2","Condition_mCherry_strict","E55")

for (k in c("Neftel_subtype","Phase","Dirks_subtype","SPollard_subtype","ccAFv2")) {
  list <- list()
  for (i in c("E17","E20","E21","E28","E31","E34","E43","E55")) {
    list[[i]] <- cherry_proportion_plot(SSE1_filtered,k,"Condition_mCherry_strict",cell_line=i)
    
  }
  p1 <- list[["E17"]]+list[["E20"]]+list[["E21"]]+list[["E28"]]+
    list[["E31"]]+list[["E34"]]+list[["E43"]]+list[["E55"]]
  ggsave(paste0("Barplot_SSE7_",k,".pdf"),p1,width = 20,height = 12)
}


# create combined figure
Combined_plotting(seurat_type="sse7",
                  plot=c("Feature_blending","Vln","UMAP"),
                  subtype=c("Condition_mCherry_strict","Neftel_subtype",
                            "ccAFv2","Phase"),
                  vln_group="Condition_mCherry_strict",
                  features=c("SOX2","SOX9"),
                  savedir="figures",
                  cell_line=F,
                  combine=T)



# statistics analysis of mcherry + cells in each subtype
percent_mCherry <- subtype_mCherry(SSE1_filtered,"Condition_mCherry_strict","all")
length(percent_mCherry)
colnames(percent_mCherry)

chisq.test(percent_mCherry[2:3,1:5]) # Neftel subtype
chisq.test(percent_mCherry[2:3,6:8]) # cell cycle Phase
chisq.test(percent_mCherry[2:3,9:10]) # Dirks subtype
chisq.test(percent_mCherry[2:3,11:12]) # SPollard subtype
chisq.test(percent_mCherry[2:3,13:20]) # cell cycle ccAFv2
chisq.test(percent_mCherry[2:3,21:28]) # cell cycle ccAF_Phase

chisq_results <- list()
for (i in unique(SSE1_filtered$cell_line)) {
  percent_mCherry <- subtype_mCherry(SSE1_filtered,"Condition_mCherry_strict",cell_line = i)
  chisq_results[[i]][["data"]] <- percent_mCherry
  chisq_results[[i]][["Neftel_subtype"]] <- chisq.test(percent_mCherry[c(1,3),1:5])
  chisq_results[[i]][["Phase"]] <- chisq.test(percent_mCherry[c(1,3),6:8])
  chisq_results[[i]][["Dirks_subtype"]] <- chisq.test(percent_mCherry[c(1,3),9:10]) 
  chisq_results[[i]][["SPollard_subtype"]] <- chisq.test(percent_mCherry[c(1,3),11:12])
  chisq_results[[i]][["ccAFv2"]] <- chisq.test(percent_mCherry[c(1,3),13:20])
}
save(chisq_results,file = "chisq_results_SSE7.rdata") # data available at https://drive.google.com/file/d/1PuyEEVdf23-46YIwIakfag-JZhMXUp4O/view?usp=drive_link


#scenic visualization
library(SCopeLoomR)
library(SCENIC)
library(VennDiagram)
library(stringr)

SSE1_filtered$line_condition <- paste(SSE1_filtered$cell_line,SSE1_filtered$Condition_mCherry_strict,sep = "_")
SSE1_filtered$line_condition <- ifelse(SSE1_filtered$cell_line %in% c("E17","E28","E43","E55"),
                                       paste0("high_",SSE1_filtered$Condition_mCherry_strict),
                                       paste0("low_",SSE1_filtered$Condition_mCherry_strict))

# function(type,compare,control="No",TSS=NULL,SSE="SSE7",rss_num=0.2,zThreshold=1,zPicking=1)
results <- SCENIC_picking("lib","Condition_mCherry_strict",TSS = "10kbp",zPicking = 0.7) # lib or line
scenic_sse7 <- SCENIC_SSE_extracting(results,"lib","Condition_mCherry_strict","SSE7")
scenic_sse7 <- scenic_sse7[which(scenic_sse7$Freq>1),]
scenic_sse7$Pvalue <- ifelse(scenic_sse7$Freq<3,0.01,scenic_sse7$Pvalue)
scenic_sse7 <- scenic_sse7[which(scenic_sse7$Pvalue<0.05),]
scenic_sse7$Pvalue <- ifelse(scenic_sse7$Freq<3,NA,scenic_sse7$Pvalue)

scenic_sse7_sub <- scenic_sse7[c("STAT1(+),","IRF9(+),","FOS(+),","MAF(+),","ETS1(+),",
                                 "JUN(+),","HIF1A(+),","SMAD1(+),","STAT2(+),","NFE2L2(+),",
                                 "SOX8(+),","SOX9(+),","SOX2(+),","BHLHE40(+),","SALL1(+),"),]

save(results,file = "SCENIC_results.RData") # data available at https://drive.google.com/file/d/1PuyEEVdf23-46YIwIakfag-JZhMXUp4O/view?usp=drive_link
## boxplot
SCENIC_Zscore_plotting(results = results,gene_df=scenic_sse7_sub,
                       cell_line = NA,
                       times=2,font_size = 8)


## Z score heatmap
expr_final <- data.frame()
expr_meta_final <- data.frame()
gene <- rownames(scenic_sse7_sub)
for (i in unique(results$results_z_score$group)) {
  expr <- results$results_z_score[which(results$results_z_score$Topic %in% gene),]
  expr <- expr[which(expr$group == i),]
  expr <- expr[,c(-4)]
  expr <- tidyr::spread(expr,Condition_mCherry_strict,Z)
  
  colnames(expr)[2:length(colnames(expr))] <- paste(i,colnames(expr)[2:length(colnames(expr))],sep = "_")
  if (length(expr_final)==0) {
    expr_final <- expr
  } else { expr_final <- merge(expr_final,expr,by="Topic",all=T)}
  
  expr_meta <- data.frame(row.names = colnames(expr)[2:length(colnames(expr))])
  expr_meta$library <- i
  expr_meta$condition <- sub(paste0(i,"_"),"",rownames(expr_meta))
  if (length(expr_meta_final)==0) {
    expr_meta_final <- expr_meta
  } else { expr_meta_final <- rbind(expr_meta_final,expr_meta)}
  
}
expr_final <- column_to_rownames(expr_final,var = "Topic")

z_score_extract <- results$results_z_score[which(results$results_z_score$Topic %in% gene),]
count <- n_distinct(z_score_extract[,2])
lib_count <- table(as.character(z_score_extract$group))
lib_count <- lib_count/count
lib_count <- as.data.frame(lib_count)

expr_meta_final_plot <- expr_meta_final[which(expr_meta_final$condition %in% c("Control_No","SSE7_Yes","SSE7_No")),]
expr_meta_final_plot$lib_order <- lib_count$Freq[match(expr_meta_final_plot$library,lib_count$Var1)]
expr_meta_final_plot <- dplyr::arrange(expr_meta_final_plot,desc(lib_order))
expr_meta_final_plot <- expr_meta_final_plot[,-3]

expr_final_plot <- expr_final[gene,rownames(expr_meta_final_plot)]

col_list <- list(library = setNames(c(brewer.pal(7,"Set2"),brewer.pal(8,"Dark2")),levels(as.factor(expr_meta_final$library))), # Phase
                 condition = c('Control_No' = "#87CEEB",
                                 'SSE7_No' = "#00008B",'SSE7_Yes' = "#ff5030"))

pheatmap::pheatmap(expr_final_plot,
                   annotation_col = expr_meta_final_plot,
                   cluster_cols = T,cluster_rows = F,
                   annotation_colors = col_list,
                   color = colorRampPalette(c("blue","white","red"))(100),
                   breaks = seq(-2, 2, length.out = 101),
                   scale = "row",na_col = "grey",
                   border_color = NA,show_colnames = F)

# scatter plot
z_score_extract <- results$results_z_score[which(results$results_z_score$Topic %in% rownames(scenic_sse7)),]

ggplot(z_score_extract,aes(y=fct_reorder(Topic, Z, .fun = max),x=Z,color=z_score_extract[,2]))+
  geom_point()+
  #geom_line(aes(group=group),position = position_dodge(width=0.3))+
  scale_color_manual(values = c('Control_No' = "#87CEEB",
                                'SSE7_No' = "#00008B",'SSE7_Yes' = "#ff5030"))+
  theme_classic() +
  theme(legend.position="none")


# order TFs by rss
a <- results$results_rss_score
a <- a[a$Topic %in% rownames(scenic_sse7),]
a <- a[a$Condition_mCherry_strict=="SSE7_Yes",]

ggplot(a) +   ggridges::geom_density_ridges(aes(x = Rss, y = Topic),alpha=0.5)

z_score_extract_wide <- tidyr::pivot_wider(a, names_from = Topic, values_from = Rss)
summary(z_score_extract_wide[,c(3:81)])

