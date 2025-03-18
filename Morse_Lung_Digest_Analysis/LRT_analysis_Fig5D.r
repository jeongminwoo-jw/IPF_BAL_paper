# Title     : LRT analysis
# Objective : LRT Test of Morse MonoMacs using DESeq2 pipeline
# Created by: JW
# Created on: 10/08/2024


library(Seurat)
library(DESeq2)
library(EnhancedVolcano)
library(paletteer)

morse <-readRDS("Morse_Myeloid_MAY2024_f.rds")
Celltype<-as.character(unique(morse$CellType_v2))

#### LRT #########

morse_sub<-morse[,!(morse$Donor %in% c("SC31DNOR","SC31NOR","SC45NOR","SC14NOR"))]
table(morse_sub$Disease_Status)

morse_sub$Disease_1<-morse_sub$Disease_Status

morse_sub$Disease_1<-gsub("IPF[(]Lower Lobe[)]","IPF_L",morse_sub$Disease_1)
morse_sub$Disease_1<-gsub("IPF_L_Mac_Deplet","IPF_L",morse_sub$Disease_1)
morse_sub$Disease_1<-gsub("IPF[(]Upper Lobe[)]","IPF_U",morse_sub$Disease_1)
morse_sub$Disease_1<-ifelse(morse_sub$Disease=="Normal","Normal",morse_sub$Disease_1)


for (i in 1:length(Celltype)) {

################ Count Matrix Generation ################

subs<-morse_sub[,morse_sub$CellType_v2 == Celltype[i]]

mtx_aggr<-AggregateExpression(subs,assays = "RNA",slot="count",group.by = "Donor")
cluster_counts<-as.data.frame(mtx_aggr$RNA)


################# Metadata generation ################

 metadata<-as.data.frame(subs@meta.data)

 meta_1<-metadata[,c("Donor","Disease_1")]

 meta_f<-unique(meta_1)

 rownames(meta_f)<-meta_f$Donor


 meta_final<-meta_f[colnames(cluster_counts),]


 ########### DESeq2 Run #########################

 meta_final$Sample_Name<-as.factor(meta_final$Donor)
 meta_final$Status<-as.factor(meta_final$Disease_1)

 dds <- DESeqDataSetFromMatrix(cluster_counts,
                               colData = meta_final,
                               design = ~ Status)


 #Prefilter
 dds <- dds[rowSums(counts(dds))> 3, ]
 rld <- rlog(dds, blind=TRUE)


 #PCA
 p1<-DESeq2::plotPCA(rld, intgroup = "Status")




  ############# LRT test ################

  dds <- DESeq(dds, test="LRT", reduced=~1)
  res_LRT <- results(dds)

  head(res_LRT[order(res_LRT$padj),], 10)

  write.table(res_LRT,file =paste0(Celltype[i],"LRT_test_result.txt",sep="\t"))



  #######Subset gene with significance
  sig_res_LRT <- res_LRT %>%
                 data.frame() %>%
                 tibble::rownames_to_column(var="gene") %>%
                 as_tibble() %>%
                 filter(padj < 0.05)

dim(sig_res_LRT)

  # Pull out sifnificant genes
  sigLRT_genes <- sig_res_LRT %>%
                  pull(gene)



  # Obtain normalized count values for those significant genes



 DEG_mat <- assay(rld)[ rownames(assay(rld))%in% sig_res_LRT$gene, colnames(assay(rld))%in% colnames(dds)]
 df <- data.frame(Group = SummarizedExperiment::colData(dds)[,c("Status")], row.names = rownames(SummarizedExperiment::colData(dds)))


 ####### Fig 5D. Heatmap_LRT test#############

 #require(lattice)
 #pdf(paste0(Celltype[i],"_DEG_LRT_Hmap_0.05.pdf"),width=8,height=6)
 #print( pheatmap::pheatmap(DEG_mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", clustering_distance_rows = "correlation",clustering_distance_cols = "correlation", fontsize_row = 1.0) )
 #dev.off()


 p<-pheatmap::pheatmap(DEG_mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", clustering_distance_rows = "correlation",clustering_distance_cols = "correlation", fontsize_row = 7.0)


  row_anno <-as.data.frame(cutree(p$tree_row, 5))
  colnames(row_anno) <- "Cluster"
  row_anno$Cluster <- as.factor(row_anno$Cluster)
  dev.off()

pal1<-paletteer_c("ggthemes::Orange-Blue Diverging", 100)

require(lattice)
pdf(paste0(Celltype[i],"_DEG_LRT_Hmap_0.05_1.pdf"),width=8,height=6)
 print(pheatmap::pheatmap(DEG_mat, cluster_rows=TRUE, color =rev(pal1),border_color = NA, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", clustering_distance_rows = "correlation",clustering_distance_cols = "correlation",annotation_row = row_anno,  fontsize_row = 1.0))
 dev.off()

 write.csv(row_anno,paste0(Celltype[i],"_Cluster_annotation.csv"))


}
