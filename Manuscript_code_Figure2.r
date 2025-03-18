library(ggplot2)
library(pheatmap)
library(Seurat)
library(SCENIC)
lilbrary(AUCell)
library(dplyr)
library(viridis)
library(scales)


################ Fig 2A. DEG Volcano Plot ##################
# (code below is repeated for each Macrophage subsets)

result<-read.delim("AM 1_DEG_0.1.txt")
result$geneid<-rownames(result)


result$FDR <- ifelse(result$padj < 0.1 , "FDR < 0.1", "ns")
result$DEG <- ifelse(result$padj < 0.1 & abs(result$log2FoldChange) > 1.0 , "FDR < 0.1 and |Log2FC| 1.0", "ns")

result1<-result[!(is.na(result$DEG)),]


ggplot(result1, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = DEG))  +ggtitle("AM1: IPF vs HC")+
        scale_color_manual(values = c("red", "grey30")) + labs(x = "Log2 Fold Change") +
        theme(legend.position = "bottom",axis.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"))+
        geom_text_repel(
              data = subset(result, DEG =="FDR < 0.1 and |Log2FC| 1.0"),
              aes(label = geneid),
              size = 3.5,
              box.padding = unit(0.35, "lines"),
              point.padding = unit(0.3, "lines")
                )
ggsave('AM1_Volcano.pdf', dpi=900,width=7, height=5, units="in")



################ Fig 2C. Heatmap (Regulon Activity score) ##################

BAL<-readRDS("BAL_obj.rds")

BAL1<-BAL[,Idents(BAL) %in% c("FABP4 AM","IGF1+ AM","CD206hi FNhi AM","monoSPP1+ AM","CXCL10+ AMs")]


regulonAUC <- importAUCfromText("aucell_result.csv")
meta<-as.data.frame(BAL1@meta.data)

meta$Group<-paste0(meta$Celltype,"_",meta$Disease)


regulonActivity_byCellType <- sapply(split(rownames(meta), meta$Group),
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



top20<-c(tail(AM1_tf[order(AM1_tf$avg_log2FC),],20)$TF,tail(AM2_tf[order(AM2_tf$avg_log2FC),],20)$TF,tail(AM3_tf[order(AM3_tf$avg_log2FC),],20)$TF,tail(CXCL10_tf[order(CXCL10_tf$avg_log2FC),],20)$TF,tail(SPP1_tf[order(SPP1_tf$avg_log2FC),],20)$TF)
union_top20<-unique(top20)
top20_final<-stringr::str_replace(union_top20, pattern="...$", "(+)")


TF_regulon<-regulonActivity_byCellType[rownames(regulonActivity_byCellType) %in% top20_final,]
TF_regulon_Scaled <- t(scale(t(TF_regulon), center = T, scale=T))



cairo_pdf("TOP_TFs_heatmap.pdf", height=10)
pheatmap::pheatmap(TF_regulon_Scaled,cluster_cols=FALSE,cluster_rows=TRUE, border_color = "black",fontsize_row=7.0)
dev.off()


################ Fig 2E. DotPlot (Regulon Activity comparison between IPF vs HC) ##################

combn1<-read.csv("TF_FCandPadj_table_Mostafavi.csv")
combn1$Regulon<-combn1$X


combn1$Regulon<-stringr::str_replace(combn1$Regulon, pattern="...$", "(+)")
combn1$X<-NULL



FC_mat<-combn1[,c(1,3,5,7,9,11)]
pval_mat<-combn1[,c(2,4,6,8,10,11)]

pcm_fc = melt(FC_mat, id = c("Regulon"))
pcm_pval = melt(pval_mat, id = c("Regulon"))



pcm_fc$variable<-gsub("_log2FC","",pcm_fc$variable)
pcm_pval$variable<-gsub("_Padj","",pcm_pval$variable)

pcm_fc$adjP<-pcm_pval$value

pcm_fc$adjP1<- -log(pcm_fc$adjP,10)
pcm_fc$adjP<-NULL

 colnames(pcm_fc)<-c("Regulon","Celltype","log2FC","Log10FDR")



  library(wesanderson)
 pal <- wes_palette("Zissou1", 100, type = "continuous")
 pal



 pcm_fc %>%
      ggplot(aes(x=Celltype, y = Regulon)) +
      geom_point( aes(color=log2FC,size= Log10FDR)) +
      theme_bw()+  labs(size = "-Log10FDR") +scale_size(range = c(1, 8))+
      theme_linedraw() + theme(panel.grid.major = element_blank()) +
      theme(axis.text.x = element_text(angle = 45, hjust=1,size=10,face="bold"),axis.text.y = element_text(size=10,face="bold"),plot.margin = margin(0.5, 0.5, 0.5, 1.0, "cm")) +
      scale_color_gradientn(colours = pal)+
      geom_vline(xintercept = seq(1.5, length(unique(pcm_fc$Regulon)) - 0.5, 1), lwd = 0.05, colour = "grey92")+
      geom_hline(yintercept = seq(1.5, length(unique(pcm_fc$Celltype)) - 0.5, 1), lwd = 0.05, colour = "grey93")


ggsave('DotPlot_TFComparison_Mostafavi.pdf', dpi=900,width=4, height=6, units="in")
