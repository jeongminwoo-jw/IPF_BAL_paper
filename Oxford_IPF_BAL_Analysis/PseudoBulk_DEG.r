# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)



#############################################################################
############### Pseudobulk analysis using DESeq2 ############################
#############################################################################

### Load merged seurat Object

BAL_obj<-readRDS("BAL_obj.RDS")

Celltype<-as.character(unique(BAL_obj$Celltype))


################ Subset the data and run DEG analysis for each Cell Types ################


for (i in 1:length(Celltype)) {


################ Count Matrix Generation ################


 clusters <- levels(as.factor(BAL_obj$Celltype))
 clusters

Idents(BAL_obj)<-BAL_obj$Celltype

table(BAL_obj$Celltype,BAL_obj$Sample)

BAL_obj_1<-subset(BAL_obj,idents=clusters[i])

Idents(BAL_obj_1)<-BAL_obj_1$Sample
BAL_small<-subset(BAL_obj_1,downsample=300)

mtx_aggr<-AggregateExpression(BAL_small,assays = "RNA",slot="count",group.by = "Sample")
cluster_counts<-as.data.frame(mtx_aggr$RNA)


################# Metadata generation ################

 metadata<-as.data.frame(BAL_small@meta.data)

 meta_1<-metadata[,c("Sample","Disease","Sex")]
 meta_f<-unique(meta_1)

 rownames(meta_f)<-meta_f$Sample


 meta_final<-meta_f[colnames(cluster_counts),]



  ########### DESeq2 Run #########################

  meta_final$Sample<-as.factor(meta_final$Sample)
  meta_final$Disease<-as.factor(meta_final$Disease)
  meta_final$Sex<-as.factor(meta_final$Sex)


  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = meta_final,
                                design = ~ 0+Status)


  #Prefilter
  dds <- dds[rowSums(counts(dds))> 3, ]
  rld <- rlog(dds, blind=TRUE)


  #PCA
  p1<-DESeq2::plotPCA(rld, intgroup = "Disease")
  p1
  p1 + theme_bw()

  ################### DEG analysis #####################


  dds$Status <- droplevels(dds$Disease)


  dds$Status  <- relevel(dds$Status , "HC" )

  dds <- DESeq(dds)
  resultsNames(dds)


  dds_status <- results(dds, contrast=c("Disease", "IPF", "HC"))
  head(dds_status)

  dds_status <- dds_status[order(abs(dds_status$stat), decreasing = TRUE),]

  dds_status <- dds_status[!grepl("^MT-", rownames(dds_status)), ]

  write.table(dds_status, file =paste0(clusters[i],"_DEG_wRiboGenes.txt",sep=""),sep="\t")

  dds_status_padj_05 <- dds_status[!is.na(dds_status$padj) & dds_status$padj < 0.05,]

  dds_status_padj_05 <- dds_status[dds_status$pvalue < 0.001,]

  ####### PseudoBulk_Heatmap
  dds_finalDEG<-dds_status_padj_05[dds_status_padj_05$log2FoldChange < -0.8 | dds_status_padj_05$log2FoldChange > 0.8,]


  DEG_mat <- assay(rld)[ rownames(assay(rld))%in% rownames(dds_finalDEG), colnames(assay(rld))%in% colnames(dds)]
  df <- data.frame(Group = SummarizedExperiment::colData(dds)[,c("Disease","Sex")], row.names = rownames(SummarizedExperiment::colData(dds)))

  require(lattice)
  pdf(paste0(clusters[i],"_DEG_Hmap_wRiboGenes.pdf",sep=""),width=10,height=8)
  print( pheatmap::pheatmap(DEG_mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", clustering_distance_rows = "correlation",clustering_distance_cols = "correlation", fontsize_row = 5.0) )
  dev.off()

   require(lattice)
   pdf(paste0(clusters[i],"_DEG_Volcano_pval_wRiboGenes.pdf",sep=""),width=9,height=7)
   print(EnhancedVolcano(dds_status,
                    lab = rownames(dds_status),
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = paste0(clusters[i]," IPF VS HC"),
                    ylab=bquote(~-Log[10]~italic(Q)),
                    pCutoff=0.05,
                    FCcutoff=0.8,labSize=4.0))
  dev.off()


  dds_status1 <- dds_status[ ! grepl('^RP[SL]', rownames(dds_status)), ]
  write.table(dds_status1, file =paste0(clusters[i],"_DEG.txt",sep=""),sep="\t")


##### Heat map generation for DEGs #############
dds_status_padj_05 <- dds_status1[!is.na(dds_status1$padj) & dds_status1$padj < 0.05,]
dds_finalDEG<-dds_status_padj_05[ abs(dds_status_padj_05$log2FoldChange) > 1.0 ,]


DEG_mat <- assay(rld)[ rownames(assay(rld))%in% rownames(dds_finalDEG), colnames(assay(rld))%in% colnames(dds)]
df <- data.frame(Group = SummarizedExperiment::colData(dds)[,c("Disease","Sex")], row.names = rownames(SummarizedExperiment::colData(dds)))

require(lattice)
pdf(paste0(clusters[i],"_DEG_Hmap.pdf",sep=""),width=10,height=8)
print( pheatmap::pheatmap(DEG_mat, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df,scale="row", clustering_distance_rows = "correlation",clustering_distance_cols = "correlation", fontsize_row = 5.0) )
dev.off()



dds_status1$FDR <- ifelse( dds_status1$padj < 0.1 , "FDR < 0.1", "ns")
dds_status1$DEG <- ifelse(result$padj < 0.1 & abs(result$log2FoldChange) > 1.0 , "FDR < 0.1 and |Log2FC| 1.0", "ns")

result1<-dds_status1[!(is.na(dds_status1$DEG)),]

############ Volcano Plot ################

require(lattice)
pdf(paste0(clusters[i],"_DEG_Hmap.pdf",sep=""),width=10,height=8)

print(  ggplot(result1, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = DEG)) + ylim(NA,7) +ggtitle(paste0(clusters[i]": IPF vs HC"))+
          scale_color_manual(values = c("red", "#2B547E")) + labs(x = "Log2 Fold Change") +
          theme(legend.position = "bottom",axis.text=element_text(size=20), plot.title = element_text(size = 15, face = "bold"))+
          geom_text_repel(
              data = subset(result, DEG =="FDR < 0.1 and |Log2FC| 1.0"),
              aes(label = geneid),
              size = 6,
              box.padding = unit(0.35, "lines"),
              point.padding = unit(0.3, "lines")
          )
 )
dev.off()


}
