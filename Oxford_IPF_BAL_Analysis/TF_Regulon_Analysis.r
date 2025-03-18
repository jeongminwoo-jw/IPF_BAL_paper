library(Seurat)
library(SCENIC)
lilbrary(AUCell)
library(dplyr)
library(viridis)
library(scales)


##########################################################################
############### TF Regulon analysis using SCENIC ##########################
###########################################################################



setwd("/ceph/project/holab/jwoo/IPF_HT/SCENIC")


#### Prepare input file to fun pyscenic
merged <- readRDS("BAL_obj.RDS")


matrix <- t( as.matrix(merged@assays$RNA@data[Matrix::rowSums(merged@assays$RNA@data) > 0, ]))
matrix <- data.frame(cell_id=rownames(matrix), matrix)
write.csv(matrix, file="matrix.all.cells.mean.cov.10.csv", row.names = F)


#### read the data back in after running pyscenic  ######
scenic <- read.csv("aucell1.csv")

rownames(scenic) <- scenic[, 1]
scenic <- t(as.matrix(scenic[, -1]))


merged[["SCENIC"]] <- CreateAssayObject(data = scenic[, Cells(merged)])

merged@active.assay <- "SCENIC"

#FeaturePlot(merged, "BATF...")
#FeaturePlot(merged, "RUNX3...")
FeaturePlot(merged, "STAT1...")+ scale_colour_viridis_c(option = "A", direction = -1)
FeaturePlot(merged, "STAT2...")+ scale_colour_viridis_c(option = "A", direction = -1)


both.cluster.reactome <- FindAllMarkers(merged,  logfc.threshold = 0, assay = "SCENIC", slot = "data",  test.use = "LR",
                                        max.cells.per.ident = 200)

write.table(both.cluster.reactome, file="TF.network.cluster.markers.xls", sep="\t", quote=FALSE, row.names = FALSE)




Idents(merged)<-merged$Celltype

for( cluster in unique(merged@active.ident)){
  print(cluster)

  mk1 <- FindMarkers(merged[, merged@active.ident == cluster], "IPF", "HC", group.by="Disease", logfc.threshold = 0, assay = "SCENIC", slot = "data",  test.use = "LR",
                                          max.cells.per.ident = 200)
  mk1$TF <- rownames(mk1)
  write.table(mk1, file=paste0(make.names(cluster), ".TF.network.cluster.diffsIPF_vs_HC.xls"), sep="\t", quote=FALSE, row.names = FALSE)

}


saveRDS(merged, file="merged.scenic.RDS")




for (tf in rownames(merged@assays$SCENIC@data) ){

  cairo_pdf(paste0(tf, "_umap_overlay.pdf"), width=8)
  print(FeaturePlot(merged, tf)+ scale_colour_viridis_c(option = "A", direction = -1) + labs(title= stringr::str_replace(tf, pattern="...$", "")) )
  dev.off()


}

for (tf in rownames(merged@assays$SCENIC@data) ){

  cairo_pdf(paste0(tf, "_umap_overlay_split.pdf"), width=11, height=5)
  print(FeaturePlot(merged, tf)+ scale_colour_viridis_c(option = "A", direction = -1) + labs(title= stringr::str_replace(tf, pattern="...$", "")) + facet_grid(.~merged$Disease) )
  dev.off()


}


for (tf in rownames(merged@assays$SCENIC@data) ){

  cairo_pdf(paste0(tf, "_violin_split.pdf"), width=10, height=10)
  print(VlnPlot(merged[, merged@active.ident %in% c("BAL Macrophages","Cycling Cells")], tf, group.by = "Disease")+ scale_colour_viridis_c(option = "A", direction = -1) + labs(title= stringr::str_replace(tf, pattern="...$", "")) + facet_wrap(~merged[, merged@active.ident %in% c("BAL Macrophages","Cycling Cells")]$Celltype) )
  dev.off()


}



############## TF_Cluster_DotPlot  ####################



TF_markers<-read.delim("TF.network.Celltype.markers.xls")

TF_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

plot <- DotPlot(merged, features=unique(top10$gene), assay = "SCENIC") + scale_colour_viridis_c(option = "D") + theme(axis.text.x = element_text(angle = 45, hjust=1))




############# Celltype Group Specific Regulons ###################

regulonAUC <- importAUCfromText("aucell1.csv")
meta<-as.data.frame(merged@meta.data)

meta_1<-meta[meta$label %in% c("BAL Macrophages","Cycling Cells"),]
meta_1$Group<-paste0(meta_1$Celltype,"_",meta_1$Disease)


regulonActivity_byCellType <- sapply(split(rownames(meta_1), meta_1$Group),
                              function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

head(regulonActivity_byCellType)




AM1_tf<-read.delim("Celltype_Differential_TF/AM_1.TF.network.cluster.diffsIPF_vs_HC.xls")
AM1_tf<-AM1_tf[AM1_tf$p_val_adj<0.05,]
AM2_tf<-read.delim("Celltype_Differential_TF/AM_2.TF.network.cluster.diffsIPF_vs_HC.xls")
AM2_tf<-AM2_tf[AM2_tf$p_val_adj<0.05,]
AM3_tf<-read.delim("Celltype_Differential_TF/AM_3.TF.network.cluster.diffsIPF_vs_HC.xls")
AM3_tf<-AM3_tf[AM3_tf$p_val_adj<0.05,]
CXCL10_tf<-read.delim("Celltype_Differential_TF/CXCL10..AM.TF.network.cluster.diffsIPF_vs_HC.xls")
CXCL10_tf<-CXCL10_tf[CXCL10_tf$p_val_adj<0.05,]
SPP1_tf<-read.delim("Celltype_Differential_TF/SPP1..AM.TF.network.cluster.diffsIPF_vs_HC.xls")
SPP1_tf<-SPP1_tf[SPP1_tf$p_val_adj<0.05,]
CyclingAM1_tf<-read.delim("Celltype_Differential_TF/Cycling.AM_1.TF.network.cluster.diffsIPF_vs_HC.xls")
CyclingAM1_tf<-CyclingAM1_tf[CyclingAM1_tf$p_val_adj<0.05,]
CyclingAM2_tf<-read.delim("Celltype_Differential_TF/Cycling.AM_2.TF.network.cluster.diffsIPF_vs_HC.xls")
CyclingAM2_tf<-CyclingAM2_tf[CyclingAM2_tf$p_val_adj<0.05,]

Naive_CD4T_tf<-read.delim("Celltype_Differential_TF/Naive.CD4..T.cells.TF.network.cluster.diffsIPF_vs_HC.xls")
Naive_CD4T_tf<-Naive_CD4T_tf[Naive_CD4T_tf$p_val_adj<0.05,]
Effector_CD4T_tf<-read.delim("Celltype_Differential_TF/Effector.CD4..T.cells.TF.network.cluster.diffsIPF_vs_HC.xls")
Effector_CD4T_tf<-Effector_CD4T_tf[Effector_CD4T_tf$p_val_adj<0.05,]
CD8T_tf<-read.delim("Celltype_Differential_TF/CD8..T.cells.TF.network.cluster.diffsIPF_vs_HC.xls")
CD8T_tf<-CD8T_tf[CD8T_tf$p_val_adj<0.05,]
NK_tf<-read.delim("Celltype_Differential_TF/NK.cells.TF.network.cluster.diffsIPF_vs_HC.xls")
NK_tf<-NK_tf[NK_tf$p_val_adj<0.05,]


top20<-c(tail(AM1_tf[order(AM1_tf$avg_log2FC),],20)$TF,tail(AM2_tf[order(AM2_tf$avg_log2FC),],20)$TF,tail(AM3_tf[order(AM3_tf$avg_log2FC),],20)$TF,tail(CXCL10_tf[order(CXCL10_tf$avg_log2FC),],20)$TF,tail(SPP1_tf[order(SPP1_tf$avg_log2FC),],20)$TF,tail(CyclingAM1_tf[order(CyclingAM1_tf$avg_log2FC),],20)$TF,tail(CyclingAM2_tf[order(CyclingAM2_tf$avg_log2FC),],20)$TF,tail(Naive_CD4T_tf[order(Naive_CD4T_tf$avg_log2FC),],20)$TF,tail(Effector_CD4T_tf[order(Effector_CD4T_tf$avg_log2FC),],20)$TF,tail(CD8T_tf[order(CD8T_tf$avg_log2FC),],20)$TF,tail(NK_tf[order(NK_tf$avg_log2FC),],20)$TF)
union_top20<-unique(top20)
top20_final<-stringr::str_replace(union_top20, pattern="...$", "(+)")

TF_regulon<-regulonActivity_byCellType[rownames(regulonActivity_byCellType) %in% top20_final,]
TF_regulon_Scaled <- t(scale(t(TF_regulon), center = T, scale=T))
clust <- pheatmap::pheatmap(TF_regulon_Scaled, color = viridisLite::viridis(option = "D", n=100), border_color = "black")



cairo_pdf("TOP_TFs_heatmap.pdf", height=10)
pheatmap::pheatmap(TF_regulon_Scaled, color = viridisLite::viridis(option = "D", n=100), border_color = "black")
dev.off()
