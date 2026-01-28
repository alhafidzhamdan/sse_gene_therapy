library(ggplot2)
library(Seurat)
library(tibble)
library(SCENIC)
library(SCopeLoomR)
library(dplyr)
library(scales)
library(stringr)
library(VennDiagram)

cherry_proportion_plot <- function(seurat_obj, mcherry, subtype, cell_line="all") {
  if (cell_line == "all") {
    a <- as.data.frame(table(seurat_obj@meta.data[,mcherry],
                             seurat_obj@meta.data[,subtype]))
  } else if (cell_line %in% c("E17","E20","E21","E28","E31","E34","E43","E55")) {
    a <- as.data.frame(table(seurat_obj@meta.data[which(seurat_obj$cell_line == cell_line),mcherry],
                             seurat_obj@meta.data[which(seurat_obj$cell_line == cell_line),subtype]))
  }
  
  color <- switch (mcherry,
                   "mCherry_count" = c('No'= '#cceeff','Yes'= '#FF88C2'),
                   "mCherry_count_strict" = c('No'= '#cceeff','Yes'= '#FF88C2'), # mCherry
                   "Neftel_subtype" = c('NPC' = '#54B345','MES' = '#ff8c8c','AC' = '#cfca50',
                                        'OPC' = '#57C3F3','hybrids' = 'grey'), # Neftel_subtype
                   "Neftel_subtype_noHybrids" = c('NPC' = '#54B345','MES' = '#ff8c8c','AC' = '#cfca50',
                                        'OPC' = '#57C3F3','hybrids' = 'grey'), # Neftel_subtype
                   "ccAF_phase" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                    'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                    'M/Early G1' = '#de425b','G1/other' = '#ea6f60','grey'), # ccAF_phase
                   "ccAFv2" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                'M/Early G1' = '#de425b','Unknown' = 'grey'), # ccAFv2
                   "Phase" = c('G1' = '#b3b261','S' = '#ffdaac','G2M' = '#f2956f'), # Phase
                   "Dirks_subtype" = c('Developmental' = "#e26565",'Injury_Response' = "#25555c"), # Dirks_subtype
                   "SPollard_subtype" = c("#339999","#cc5d8a") # SPollard_subtype
  )
  
  a$proportion <- NA
  
  for (i in unique(a$Var2)) {
    b <- a[which(a$Var2 %in% i),]
    b$Var1 <- b$Freq/sum(b$Freq)
    a[which(a$Var2 %in% i),4] <- b$Var1
  }
  
  this_tile <- paste0(cell_line,
                      "\nComposition of ",subtype," in ",mcherry)
  
  if (mcherry %in% c("ccAFv2","ccAF_phase")) {
    a$Var1 <- factor(a$Var1,levels = c("Neural G0","G1","Late G1","S","S/G2","G2/M","M/Early G1","G1/other","Unknown"))
  } else if (mcherry=="Phase") {
    a$Var1 <- factor(a$Var1,levels = c("G1","S","G2M"))
  } else if (mcherry=="Neftel_subtype") {
    a$Var1 <- factor(a$Var1,levels = c("MES","AC","NPC","OPC","hybrids"))
  }
  
  p1 <- ggplot(a,aes(x=Var2,y=proportion, fill = Var1))+
    geom_col(position="stack",width = 0.5)+ # stack:堆叠图
    scale_y_continuous(expand=c(0,0))+ # 调整y轴属性，使柱子与x轴坐标接触
    scale_fill_manual(values = color)+
    labs(x='Samples',y= paste0(subtype,'(%)'))+
    ggtitle(this_tile)+
    theme(plot.title = element_text(size=15,hjust = 0.5))+ # 设置X轴和Y轴信息
    theme(panel.grid = element_blank(),panel.background = element_rect(fill='transparent'),strip.text = element_text(size=12))+ 
    theme(axis.text.x = element_text(angle = 45, hjust = .5, vjust = .5),axis.text = element_text(size = 12), axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))
  
  return(p1)
}

subtype_plot <- function(seurat_obj, subtype, cell_line="all", condition="Control_No") {
  print(condition)
  if (condition %in% unique(seurat_obj$Condition_mCherry_strict)) {
    seurat_obj <- subset(seurat_obj,Condition_mCherry_strict==condition)
  }
  
  if (cell_line == "all") {
    a <- as.data.frame(table(seurat_obj@meta.data[,subtype],
                             seurat_obj@meta.data[,"cell_line"]))
  } else if (cell_line %in% unique(seurat_obj$cell_line)) {
    a <- as.data.frame(table(seurat_obj@meta.data[which(seurat_obj$cell_line == cell_line),subtype],
                             seurat_obj@meta.data[which(seurat_obj$cell_line == cell_line),"cell_line"]))
  }
  
  color <- switch (subtype,
                   "mCherry_count" = c('#cceeff','#FF88C2'),
                   "mCherry_count_strict" = c('#cceeff','#FF88C2'), # mCherry
                   "Neftel_subtype" = c('NPC' = '#54B345','MES' = '#ff8c8c','AC' = '#cfca50',
                                        'OPC' = '#57C3F3','hybrids' = 'grey'), # Neftel_subtype
                   "Neftel_subtype_noHybrids" = c('NPC' = '#54B345','MES' = '#ff8c8c','AC' = '#cfca50',
                                        'OPC' = '#57C3F3','hybrids' = 'grey'), # Neftel_subtype
                   "ccAF_phase" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                    'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                    'M/Early G1' = '#de425b','G1/other' = '#ea6f60','grey'), # ccAF_phase
                   "ccAFv2" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                'M/Early G1' = '#de425b','Unknown' = 'grey'), # ccAFv2
                   "Phase" = c('G1' = '#b3b261','S' = '#ffdaac','G2M' = '#f2956f'), # Phase
                   "Dirks_subtype" = c('Developmental' = "#e26565",'Injury_Response' = "#25555c"), # Dirks_subtype
                   "SPollard_subtype" = c("#339999","#cc5d8a") # SPollard_subtype
  )
  
  this_tile <- paste0(cell_line,
                      '\nThe number of cells is ',sum(a$Freq))
  
  if (subtype %in% c("ccAFv2","ccAF_phase")) {
    a$Var1 <- factor(a$Var1,levels = c("Neural G0","G1","Late G1","S","S/G2","G2/M","M/Early G1","G1/other","Unknown"))
  } else if (subtype=="Phase") {
    a$Var1 <- factor(a$Var1,levels = c("G1","S","G2M"))
  } else if (subtype=="Neftel_subtype") {
    a$Var1 <- factor(a$Var1,levels = c("MES","AC","NPC","OPC","hybrids"))
  }
  
  if (cell_line == "all") {
    
    a$proportion <- NA
    for (i in unique(a$Var2)) {
      b <- a[which(a$Var2 %in% i),]
      b$Var1 <- b$Freq/sum(b$Freq)
      a[which(a$Var2 %in% i),4] <- b$Var1
    }
    
    p1 <- ggplot(a,aes(x=Var2,y=round(proportion*100,1), fill = Var1))+
      geom_col(position="stack",width = 0.5)+ # stack:堆叠图
      scale_y_continuous(expand=c(0,0))+ # 调整y轴属性，使柱子与x轴坐标接触
      scale_fill_manual(values = color)+
      labs(x='Cell lines',y=paste0(subtype,'(%)'))+
      theme(plot.title = element_text(size=15,hjust = 0.5))+ # 设置X轴和Y轴信息
      theme(panel.grid = element_blank(),panel.background = element_rect(fill='transparent'),strip.text = element_text(size=12))+ 
      theme(axis.text.x = element_text(angle = 45, hjust = .5, vjust = .5),axis.text = element_text(size = 12), 
            axis.title = element_text(size = 13), legend.title = element_blank(), legend.text = element_text(size = 11))
  } else if (cell_line %in% c("E17","E20","E21","E28","E31","E34","E43","E55")) {
    
    a[,"ratio"] <- percent(a$Freq/sum(a$Freq),0.01)
    
    p1 <- ggplot(a, aes(x="", y=Freq, fill=a[,1]))+
      geom_bar(width = 1, stat = "identity")+coord_polar(theta = 'y')+
      labs(x = '', y = '', title = '')+
      theme(axis.text = element_blank())+
      theme(axis.ticks = element_blank())+theme_void()+
      geom_text_repel(aes(x = 1.2, label = a[,"ratio"]),
                      position = position_stack(vjust = 0.5), color = "black")+
      ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
      scale_fill_manual(values = color)
  }
  return(p1)
}

Combined_plotting <- function(seurat_type,plot,subtype,features,vln_group,
                              savedir,cell_line=F,combine=T){
  
  path <- paste0(getwd(),"/",savedir)
  if (!dir.exists(path)) {
    dir.create(path)
  }
  print(paste0("figures are save in ",path))
  
  plot_umap <- function(seurat_object, subtype) {
    umap <- list()
    for (pheno in subtype) {
      cols <- switch (pheno,
                      "mCherry_count" = c('No'= '#cceeff','Yes'= '#FF88C2'),
                      "mCherry_count_strict" = c('No'= '#cceeff','Yes'= '#FF88C2'),
                      "Neftel_subtype" = c('NPC' = '#54B345','MES' = '#ff8c8c','AC' = '#cfca50',
                                           'OPC' = '#57C3F3','hybrids' = 'grey'),
                      "ccAF_phase" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                       'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                       'M/Early G1' = '#de425b','G1/other' = '#ea6f60','grey'),
                      "ccAFv2" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                   'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                   'M/Early G1' = '#de425b','Unknown' = 'grey'), 
                      "Phase" = c('G1' = '#b3b261','S' = '#ffdaac','G2M' = '#f2956f'), 
                      "Dirks_subtype" = c('Developmental' = "#e26565",'Injury_Response' = "#25555c"),
                      "SPollard_subtype" = c("#339999","#cc5d8a"),
                      "Condition_mCherry_strict" = c('Control_No' = "#87CEEB",'SSE7_No' = "#00008B",'SSE7_Yes' = "#ff5030")
      )
      
      umap[[pheno]] <- DimPlot(seurat_object,
                               reduction = "umap",#umap和tsne两种展示方法选一种
                               group.by = pheno,#自己定义展示的分组类别
                               pt.size = 0.1,cols = cols,label = F,
                               order = c("SSE7_Yes","SSE7_No","Control_No"))
    }
    return(umap)
  }
  
  seurat_path <- switch(seurat_type,
                        "org" = paste0(".rds"),
                        "ctr" = paste0("_ctr.rds"),
                        "sse7" = paste0("_sse7.rds"),
                        "sse9" = paste0("_sse9.rds")
  )
  
  if (!cell_line) {
    seurat_object <- readRDS(paste0("SSE1_5_seurat_noR_annotation",seurat_path))
    
    result_list <- list()
    if ("Feature_blending" %in% plot) {
      p <- FeaturePlot(seurat_object, features,
                       col=c("lightgrey", "#ff0000", "#008500"),
                       blend = T,blend.threshold = 0,order = T,
                       pt.size = 1)
      result_list[["Feature_blending"]] <- p
      ggsave(paste0(path,"/Feature_blending_all_",features[1],"_",features[2],".pdf"),p,width = 14,height = 5)
      print("Feature_blending plot is done")
    }
    
    if ("Vln" %in% plot) {
      color <- switch (vln_group,
                       "mCherry_count" = c('#cceeff','#FF88C2'),
                       "mCherry_count_strict" = c('#cceeff','#FF88C2'), # mCherry
                       "Condition_mCherry_strict" = c('Control_No' = "#87CEEB",'SSE7_No' = "#00008B",'SSE7_Yes' = "#ff5030"), # mCherry
                       "Neftel_subtype" = c('NPC' = '#54B345','MES' = '#ff8c8c','AC' = '#cfca50',
                                            'OPC' = '#57C3F3','hybrids' = 'grey'), # Neftel_subtype
                       "ccAF_phase" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                        'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                        'M/Early G1' = '#de425b','G1/other' = '#ea6f60','grey'), # ccAF_phase
                       "ccAFv2" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                    'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                    'M/Early G1' = '#de425b','Unknown' = 'grey'), # ccAFv2
                       "Phase" = c('G1' = '#b3b261','S' = '#ffdaac','G2M' = '#f2956f'), # Phase
                       "Dirks_subtype" = c('Developmental' = "#e26565",'Injury_Response' = "#25555c"), # Dirks_subtype
                       "SPollard_subtype" = c("#339999","#cc5d8a") # SPollard_subtype
      )
      p <- VlnPlot(seurat_object, features,pt.size = 0,group.by = vln_group,cols = color)
      result_list[["Vln"]] <- p
      ggsave(paste0(path,"/Vln_all_",features[1],"_",features[2],".pdf"),p,width = 14)
      print("Vln plot is done")
    }
    
    if ("UMAP" %in% plot) {
      result_list[["UMAP"]]  <- plot_umap(seurat_object, subtype)
      for (umapnames in subtype) {
        ggsave(paste0(path,"/UMAP_",line,"_",umapnames,".pdf"),result_list[['UMAP']][[umapnames]],height = 5)
      }
      print("UMAP plot is done")
    }
    
  } else{
  
  combined_list <- list()
  for (line in c("E17","E20","E21","E28","E31","E34","E43","E55")) {
    seurat_object <- readRDS(paste0("SSE1_6_seurat_",line,seurat_path))
    print(paste0(line," is processing"))
    
    if ("Feature_blending" %in% plot) {
      p <- FeaturePlot(seurat_object, features,
                       col=c("#cfdfd9","#e66276","#834eae"),
                       blend = T,blend.threshold = 0,
                       pt.size = 0.02)
      combined_list[[line]][["Feature_blending"]] <- p
      ggsave(paste0(path,"/Feature_blending_",line,"_",features[1],"_",features[2],".pdf"),p,width = 14,height = 5)
      print("Feature_blending plot is done")
    }
    
    if ("Vln" %in% plot) {
      color <- switch (vln_group,
                       "mCherry_count" = c('#cceeff','#FF88C2'),
                       "mCherry_count_strict" = c('#cceeff','#FF88C2'), # mCherry
                       "Condition_mCherry_strict" = c('Control_No' = "#87CEEB",'SSE7_No' = "#00008B",'SSE7_Yes' = "#ff5030"), # mCherry
                       "Neftel_subtype" = c('NPC' = '#54B345','MES' = '#ff8c8c','AC' = '#cfca50',
                                            'OPC' = '#57C3F3','hybrids' = 'grey'), # Neftel_subtype
                       "ccAF_phase" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                        'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                        'M/Early G1' = '#de425b','G1/other' = '#ea6f60','grey'), # ccAF_phase
                       "ccAFv2" = c('Neural G0' = '#488f31', 'G1' = '#93a74d', 'Late G1' = '#cfbf78',
                                    'S' = '#fedaac', 'S/G2' = '#f2ab7d', 'G2/M' = '#e5795e', 
                                    'M/Early G1' = '#de425b','Unknown' = 'grey'), # ccAFv2
                       "Phase" = c('G1' = '#b3b261','S' = '#ffdaac','G2M' = '#f2956f'), # Phase
                       "Dirks_subtype" = c('Developmental' = "#e26565",'Injury_Response' = "#25555c"), # Dirks_subtype
                       "SPollard_subtype" = c("#339999","#cc5d8a") # SPollard_subtype
      )
      p <- VlnPlot(seurat_object, features,pt.size = 0,group.by = vln_group,cols = color)
      combined_list[[line]][["Vln"]] <- p
      ggsave(paste0(path,"/Vln_",line,"_",features[1],"_",features[2],".pdf"),p,width = 14)
      print("Vln plot is done")
    }
    
    if ("UMAP" %in% plot) {
      combined_list[[line]][["UMAP"]]  <- plot_umap(seurat_object, subtype)
      for (umapnames in subtype) {
        ggsave(paste0(path,"/UMAP_",line,"_",umapnames,".pdf"),combined_list[[line]][['UMAP']][[umapnames]],height = 5)
      }
      print("UMAP plot is done")
    }
  }
  if(combine){
    print("UMAP and Vln are combining")
    
    if ("UMAP" %in% plot) {
      for (sub_umap in subtype) {
        p <- combined_list[["E17"]][['UMAP']][[sub_umap]]+combined_list[["E20"]][['UMAP']][[sub_umap]]+combined_list[["E21"]][['UMAP']][[sub_umap]]+
          combined_list[["E28"]][['UMAP']][[sub_umap]]+combined_list[["E31"]][['UMAP']][[sub_umap]]+combined_list[["E34"]][['UMAP']][[sub_umap]]+
          combined_list[["E43"]][['UMAP']][[sub_umap]]+combined_list[["E55"]][['UMAP']][[sub_umap]]
        ggsave(paste0(path,"/UMAP_combined_",sub_umap,".pdf"),p,height = 12,width = 20)}
    }
    
    if ("Vln" %in% plot) {
      p <- combined_list[["E17"]][["Vln"]]+combined_list[["E20"]][["Vln"]]+combined_list[["E21"]][["Vln"]]+
        combined_list[["E28"]][["Vln"]]+combined_list[["E31"]][["Vln"]]+combined_list[["E34"]][["Vln"]]+
        combined_list[["E43"]][["Vln"]]+combined_list[["E55"]][["Vln"]]
      ggsave(paste0(path,"/Vln_combined.pdf"),p,height = 12,width = 20)
    }
  }
  }
}  


SCENIC_picking <- function(type,compare,TSS=NULL,SSE="SSE7",zPicking=1){
  print(paste0("TFs selection: z score > ",zPicking))
  
  path <- switch (type,
                  "lib" = "./scenic/library",
                  "line" = "./scenic/line"
  )
  
  if (!is.null(TSS)) {
    path <- switch (SSE,
                    "SSE7" = paste0(path,"_sse7"),
                    "SSE9" = paste0(path,"_sse9"))
    
    path <- switch (TSS,
                    "10kbp" = paste0(path,"_10kbp/"),
                    "500bp" = paste0(path,"_500bp/"))
  } 
  
  print(paste0("Comparison group is: ",compare))
  
  group <- switch (type,
                   "lib" = paste0("lib", 1:15),
                   "line" = unique(SSE1_filtered$cell_line)
  )
  
  data = "_regulon_AUC.loom"
  
  print(paste0("data path is: ",path))
  
  SCENIC_picking_list <- list()
  condition_strict_z <- list()
  regulon_list <- list()
  rss_list <- list()
  rss_final <- data.frame()
  z_score_final <- data.frame()
  threshold_list <- list()
  AUC_list <- list()
  for (i in group) {
    AUC_lib <- open_loom(paste0(path,i,data))
    
    regulon_incidMat <- get_regulons(AUC_lib,column.attr.name = "Regulons")
    regulons <- regulonsToGeneLists(regulon_incidMat)
    regulon_list[[i]] <- names(regulons)
    
    regulon_AUC <- get_regulons_AUC(AUC_lib,column.attr.name = "RegulonsAUC")
    AUC_list[[i]] <- regulon_AUC
    regulonAucThresholds <- get_regulon_thresholds(AUC_lib)
    threshold_change <- setNames(names(regulonAucThresholds),regulonAucThresholds)
    threshold_list[[i]] <- threshold_change
    
    close_loom(AUC_lib)
    
    regulon_AUC@metadata <- SSE1_filtered@meta.data[which(rownames(SSE1_filtered@meta.data) %in% colnames(regulon_AUC)),]
    
    rss <- calcRSS(AUC = regulon_AUC,
                   cellAnnotation = regulon_AUC@metadata[,compare])
    rss_list[[i]] <- rss
    rss <- na.omit(rss)
    
    rss_temp <- reshape2::melt(rss)
    colnames(rss_temp) <- c("Topic","Condition_mCherry_strict","Rss")
    rss_temp$group <- i
    
    if (length(rss_final)==0) {
      rss_final <- rss_temp
    } else {rss_final <- rbind(rss_final,rss_temp)}
    
    a <- rss  %>% scale %>% as.data.frame %>% rownames_to_column(var = "Topic")
    
    suppressMessages({
    a <- reshape2::melt(a)
    })
    
    colnames(a)[2:3] <- c("Condition_mCherry_strict","Z")
    a$group <- i
    
    for (j in 1:length(colnames(rss))) {
      b <- a[which(a[,2]==colnames(rss)[j] & a$Z>zPicking),]
      condition_strict_z[[colnames(rss)[j]]][[i]] <- as.character(b$Topic)
    }
    
    if (length(z_score_final)==0) {
      z_score_final <- a
    } else {z_score_final <- rbind(z_score_final,a)}
  }
  
  SCENIC_picking_list[['results_z_list']] <- condition_strict_z
  SCENIC_picking_list[['results_regulon_names']] <- regulon_list
  SCENIC_picking_list[['results_raw_rss']] <- rss_list
  SCENIC_picking_list[['results_AUC']] <- AUC_list
  SCENIC_picking_list[['results_threshold']] <- threshold_list
  SCENIC_picking_list[['results_z_score']] <- z_score_final
  SCENIC_picking_list[['results_rss_score']] <- rss_final
  
  return(SCENIC_picking_list)
}

SCENIC_SSE_extracting <- function(results=SCENIC_picking,type,compare,sse){
  
  if (compare=="Condition_mCherry_strict") {
    group <- NA
  } else {
    group <- unique(sub("_.*", "", SSE1_filtered@meta.data[,compare]))
  }
  
  result_list <- list()
  if (any(is.na(group))) {
    ctr <- "Control_No"
    sse_n <- paste(sse,"No",sep="_")
    sse_y <- paste(sse,"Yes",sep="_")
    
    if (type=="lib") {
      
      scenic_ctr <- table(unlist(results$results_z_list[[ctr]])) %>% names
      scenic_sse_n <- table(unlist(results$results_z_list[[sse_n]])) %>% names
      scenic_sse_y <- table(unlist(results$results_z_list[[sse_y]])) %>% names
      
      if (is.null(scenic_sse_y)) {
        print(paste0(sse_y,": no TFs"))
        return("adjust zPicking in SCENIC_picking")
      } else { print(paste0(sse_y,": ",length(scenic_sse_y)," TFs"))}
      
      if (is.null(scenic_ctr)) {
        print("Control: no TFs")
      } else { print(paste0("Control: ",length(scenic_ctr)," TFs"))}
      
      if (is.null(scenic_sse_n)) {
        print(paste0(sse_n,": no TFs"))
      } else { print(paste0(sse_n,": ",length(scenic_sse_n)," TFs"))}
      
      z_score_extract <- results$results_z_score[which(results$results_z_score$Topic %in% scenic_sse_y),]
      count <- n_distinct(z_score_extract[,2])
      TF_count <- table(as.character(z_score_extract$Topic))
      TF_count <- TF_count/count
      TF_count <- as.data.frame(TF_count)
      TF_count <- dplyr::arrange(TF_count,desc(Freq))
      TF_count$Pvalue <- NA
      TF_count$sse_n_vs_ctr <- NA
      TF_count$sse_y_vs_ctr <- NA
      TF_count$sse_y_vs_sse_n <- NA
      TF_count$control <- NA
      TF_count$sse_n <- NA
      TF_count$sse_y <- NA
      TF_count <- column_to_rownames(TF_count,var = "Var1")
      for (TF in unique(as.character(z_score_extract$Topic))) {
        anova_df <- z_score_extract[which(z_score_extract$Topic == TF),]
        TF_count[TF,"control"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict==ctr)])
        TF_count[TF,"sse_n"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict==sse_n)])
        TF_count[TF,"sse_y"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict==sse_y)])
        
        if (TF %in% rownames(TF_count)[which(TF_count$Freq<3)]) {
          TF_count[TF,"Pvalue"] <- 0
          next
        }
        
        anova_results <- aov(Z ~ Condition_mCherry_strict + Error(group/Condition_mCherry_strict),anova_df)
        anova_summary <- summary(anova_results)
        
        fit <- lm(Z ~ Condition_mCherry_strict + group, data = anova_df)
        emm <- emmeans(fit, ~ Condition_mCherry_strict)
        emm_results <- pairs(emm, adjust = "bonferroni") %>% as.data.frame
        emm_results <- column_to_rownames(emm_results,var = "contrast")
        
        TF_count[TF,"Pvalue"] <- anova_summary[[2]][[1]][["Pr(>F)"]][1]
        TF_count[TF,"sse_n_vs_ctr"] <- emm_results[paste(ctr,sse_n,sep = " - "),"p.value"]
        TF_count[TF,"sse_y_vs_ctr"] <- emm_results[paste(ctr,sse_y,sep = " - "),"p.value"]
        TF_count[TF,"sse_y_vs_sse_n"] <- emm_results[paste(sse_n,sse_y,sep = " - "),"p.value"]
      }
      TF_count <- TF_count[which(TF_count$sse_y > TF_count$control &
                                   TF_count$sse_y > TF_count$sse_n),]
      
      TF_count[rownames(TF_count)[which(TF_count$Freq<3)],"Pvalue"] <- NA
      
      return(TF_count)
      
    } else {
      cell_line_extracting <- list()
      for (i in unique(SSE1_filtered$cell_line)) {
        print(i)
        
        scenic_ctr <- table(unlist(results[[i]]$results_z_list[[ctr]])) %>% names
        scenic_sse_n <- table(unlist(results[[i]]$results_z_list[[sse_n]])) %>% names
        scenic_sse_y <- table(unlist(results[[i]]$results_z_list[[sse_y]])) %>% names
        
        if (is.null(scenic_sse_y)) {
          print(paste0(sse_y,": no TFs"))
          return("adjust zPicking in SCENIC_picking")
        } else { print(paste0(sse_y,": ",length(scenic_sse_y)," TFs"))}
        
        if (is.null(scenic_ctr)) {
          print("Control: no TFs")
        } else { print(paste0("Control: ",length(scenic_ctr)," TFs"))}
        
        if (is.null(scenic_sse_n)) {
          print(paste0(sse_n,": no TFs"))
        } else { print(paste0(sse_n,": ",length(scenic_sse_n)," TFs"))}
        
        z_score_extract <- results[[i]]$results_z_score[which(results[[i]]$results_z_score$Topic %in% scenic_sse_y),]
        count <- n_distinct(z_score_extract[,2])
        TF_count <- table(as.character(z_score_extract$Topic))
        TF_count <- TF_count/count
        TF_count <- as.data.frame(TF_count)
        TF_count <- dplyr::arrange(TF_count,desc(Freq))
        TF_count$Pvalue <- NA
        TF_count$sse7n_vs_ctr <- NA
        TF_count$sse7y_vs_ctr <- NA
        TF_count$sse7y_vs_sse7n <- NA
        TF_count$control <- NA
        TF_count$sse7_n <- NA
        TF_count$sse7_y <- NA
        TF_count <- column_to_rownames(TF_count,var = "Var1")
        for (TF in unique(as.character(z_score_extract$Topic))) {
          anova_df <- z_score_extract[which(z_score_extract$Topic == TF),]
          anova_df$Condition_mCherry_strict <- factor(
            anova_df$Condition_mCherry_strict,
            levels = c("Control_No", "SSE7_No", "SSE7_Yes")   # fixed order
          )
          TF_count[TF,"control"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict=="Control_No")])
          TF_count[TF,"sse7_n"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict=="SSE7_No")])
          TF_count[TF,"sse7_y"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict=="SSE7_Yes")])
          
          if (TF %in% rownames(TF_count)[which(TF_count$Freq<3)]) {
            TF_count[TF,"Pvalue"] <- 0
            next
          }
          
          anova_results <- aov(Z ~ Condition_mCherry_strict + Error(group/Condition_mCherry_strict),anova_df)
          anova_summary <- summary(anova_results)
          
          fit <- lm(Z ~ Condition_mCherry_strict + group, data = anova_df)
          emm <- emmeans(fit, ~ Condition_mCherry_strict)
          emm_results <- pairs(emm, adjust = "bonferroni") %>% as.data.frame
          emm_results <- column_to_rownames(emm_results,var = "contrast")
          
          TF_count[TF,"Pvalue"] <- anova_summary[[2]][[1]][["Pr(>F)"]][1]
          TF_count[TF,"sse7n_vs_ctr"] <- emm_results["Control_No - SSE7_No","p.value"]
          TF_count[TF,"sse7y_vs_ctr"] <- emm_results["Control_No - SSE7_Yes","p.value"]
          TF_count[TF,"sse7y_vs_sse7n"] <- emm_results["SSE7_No - SSE7_Yes","p.value"]
        }
        TF_count <- TF_count[which(TF_count$sse7_y > TF_count$control &
                                     TF_count$sse7_y > TF_count$sse7_n &
                                     TF_count$sse7_y > 0.7),]
        
        TF_count[rownames(TF_count)[which(TF_count$Freq<3)],"Pvalue"] <- NA
        
        cell_line_extracting[[i]] <- TF_count
      }
      
      print(paste0("data focuse on: sse7 in each line"))
      return(cell_line_extracting)
    }
    
  } else {
    for (i in group) {
      ctr <- paste(i,"Control_No",sep="_")
      sse_n <- paste(i,sse,"No",sep="_")
      sse_y <- paste(i,sse,"Yes",sep="_")
      
      scenic_ctr <- table(unlist(results$results_z_list[[ctr]])) %>% names
      scenic_sse_n <- table(unlist(results$results_z_list[[sse_n]])) %>% names
      scenic_sse_y <- table(unlist(results$results_z_list[[sse_y]])) %>% names
      
      
      if (is.null(scenic_sse_y)) {
        print(paste0(sse_y,": no TFs"))
        return("adjust zPicking in SCENIC_picking")
      } else { print(paste0(sse_y,": ",length(scenic_sse_y)," TFs"))}
      
      if (is.null(scenic_ctr)) {
        print("Control: no TFs")
      } else { print(paste0("Control: ",length(scenic_ctr)," TFs"))}
      
      if (is.null(scenic_sse_n)) {
        print(paste0(sse_n,": no TFs"))
      } else { print(paste0(sse_n,": ",length(scenic_sse_n)," TFs"))}
      
      
      z_score_extract <- results$results_z_score[which(results$results_z_score$Topic %in% scenic_sse_y),]
      colnames(z_score_extract)[2] <- "Condition_mCherry_strict"
      z_score_extract <- z_score_extract[which(sub("_.*", "", z_score_extract$Condition_mCherry_strict) == i),]
      z_score_extract$Condition_mCherry_strict <- sub("^[^_]+_", "", z_score_extract$Condition_mCherry_strict)
      
      count <- n_distinct(z_score_extract[,2])
      TF_count <- table(as.character(z_score_extract$Topic))
      TF_count <- TF_count/count
      TF_count <- as.data.frame(TF_count)
      TF_count <- dplyr::arrange(TF_count,desc(Freq))
      TF_count$Pvalue <- NA
      TF_count$sse7n_vs_ctr <- NA
      TF_count$sse7y_vs_ctr <- NA
      TF_count$sse7y_vs_sse7n <- NA
      TF_count$control <- NA
      TF_count$sse7_n <- NA
      TF_count$sse7_y <- NA
      TF_count <- column_to_rownames(TF_count,var = "Var1")
      for (TF in unique(as.character(z_score_extract$Topic))) {
        anova_df <- z_score_extract[which(z_score_extract$Topic == TF),]
        TF_count[TF,"control"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict=="Control_No")])
        TF_count[TF,"sse7_n"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict=="SSE7_No")])
        TF_count[TF,"sse7_y"] <- mean(anova_df$Z[which(anova_df$Condition_mCherry_strict=="SSE7_Yes")])
        
        if (TF %in% rownames(TF_count)[which(TF_count$Freq<3)]) {
          TF_count[TF,"Pvalue"] <- 0
          next
        }
        
        anova_results <- aov(Z ~ Condition_mCherry_strict + Error(group/Condition_mCherry_strict),anova_df)
        anova_summary <- summary(anova_results)
        
        fit <- lm(Z ~ Condition_mCherry_strict + group, data = anova_df)
        emm <- emmeans(fit, ~ Condition_mCherry_strict)
        emm_results <- pairs(emm, adjust = "bonferroni") %>% as.data.frame
        emm_results <- column_to_rownames(emm_results,var = "contrast")
        
        TF_count[TF,"Pvalue"] <- anova_summary[[2]][[1]][["Pr(>F)"]][1]
        TF_count[TF,"sse7n_vs_ctr"] <- emm_results["Control_No - SSE7_No","p.value"]
        TF_count[TF,"sse7y_vs_ctr"] <- emm_results["Control_No - SSE7_Yes","p.value"]
        TF_count[TF,"sse7y_vs_sse7n"] <- emm_results["SSE7_No - SSE7_Yes","p.value"]
      }
      TF_count <- TF_count[which(TF_count$sse7_y > TF_count$control &
                                   TF_count$sse7_y > TF_count$sse7_n &
                                   TF_count$sse7_y > 0.7),]
      
      TF_count[rownames(TF_count)[which(TF_count$Freq<3)],"Pvalue"] <- NA
      
      result_list[[i]] <- TF_count
    }
    return(result_list)
  }
}

SCENIC_Zscore_plotting <- function(SCENIC_picking,results=results,gene_df,times,
                                   font_size=8,cell_line=NA){
  gene <- rownames(gene_df)
  
  if (!is.na(cell_line)) {
    z_score_extract <- results$results_z_score[which(substr(results$results_z_score$line_condition,1,3) == cell_line),]
    z_score_extract <- z_score_extract[which(z_score_extract$Topic %in% gene),]
  } else {
    z_score_extract <- results$results_z_score[which(results$results_z_score$Topic %in% gene),]
  }
    
  
  print("TFs' distribution in libraries: ")
  print(table(gene_df[,1]))
  
  group_count_name <- rownames(gene_df)[which(gene_df$Freq>=times)]
  z_score_extract <- z_score_extract[which(z_score_extract$Topic %in% group_count_name),]
  if (times>0) {
    print(paste0("Display TFs in more than ", times," libraries"))
  } else {print("Display all TFs")}
  

  p1 <- ggplot(z_score_extract,aes(x=fct_reorder(Topic, Z, .fun = max),y=Z,color=z_score_extract[,2]))+
    geom_jitter(alpha=0.8,
                position=position_jitterdodge(jitter.width = 0.3, 
                                              jitter.height = 0, 
                                              dodge.width = 0.8))+
    geom_boxplot(alpha=0.2,width=0.3,
                 position=position_dodge(width=0.8),
                 size=0.7,outlier.colour = NA)+
    scale_color_manual(values = c('Control_No' = "#87CEEB",
                                  'SSE7_No' = "#00008B",'SSE7_Yes' = "#ff5030"))+
    theme_classic() +
    theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle=45,size = font_size,hjust = 1,vjust = 1))
  
  return(p1)
}

binarizeAUC <- function(auc, thresholds){
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}


