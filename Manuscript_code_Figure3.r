require(Seurat)
require(dplyr)
require(Matrix)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(Hmisc)
library(RColorBrewer)
library(SCENIC)
lilbrary(AUCell)


######## Fig 3A. UMAP (Oxford PBMC)#########

Blood<-readRDS("Blood_obj.rds")

print(head(Blood@meta.data))
DefaultAssay(object = Blood) <- "RNA"
Idents(Blood)<-Blood@meta.data$CellType
DimPlot(Blood, reduction = "umap", group.by="CellType",label=T,repel = T)
ggsave('IPF_Blood_UMAP.pdf', dpi=900,width=7, height=5, units="in")


######## Fig 3B. UMAP (Monocyte Re-clustering )#########

Mono<-Blood[,Blood$CellType %in% c("cMono","ncMono")]



###prelim clustering
Mono@active.assay <- "RNA"
Mono <- NormalizeData(Mono)
Mono <- FindVariableFeatures(Mono)
Mono <- ScaleData(Mono)
Mono <- RunPCA(object = Mono, verbose = FALSE)
ElbowPlot(Mono, ndims = 50)
Mono <- RunUMAP(object = Mono, dims = 1:20, verbose = FALSE)
Mono <- FindNeighbors(object = Mono, dims = 1:20, verbose = FALSE)
Mono <- FindClusters(object = Mono, verbose = FALSE, resolution = .7)
DimPlot(object = Mono, label = TRUE, reduction="umap")

DimPlot(object = Mono, group.by="orig.ident", reduction="umap")


Mono <- harmony::RunHarmony(Mono, group.by.vars = "orig.ident")
Mono <- RunUMAP(object = Mono, dims = 1:20, verbose = FALSE, reduction = "harmony",reduction.name = "harmony.umap")

Mono <- FindNeighbors(object = Mono, dims = 1:20, verbose = FALSE,reduction = "harmony",reduction.name = "harmony.umap")
Mono <- FindClusters(object = Mono, verbose = FALSE, resolution =0.7)

 DimPlot(object = Mono, label = TRUE, reduction="harmony.umap")

Idents(Mono)<-Mono$RNA_snn_res.0.7
markers <- FindAllMarkers(Mono, max.cells.per.ident = 200)


Mono <- RenameIdents(Mono, "0"="cMono S100A8/9/12Hi")
Mono <- RenameIdents(Mono, "1"="cMono S100A8/9/12Hi")
Mono <- RenameIdents(Mono, "2"="cMono MHCIIhi")
Mono <- RenameIdents(Mono, "3"="NcMono")
Mono <- RenameIdents(Mono, "4"="Mono-T Cell Doublets")
Mono <- RenameIdents(Mono, "6"="Mono-T Cell Doublets")
Mono <- RenameIdents(Mono, "5"="Type I IFN Mono")
Mono <- RenameIdents(Mono, "7"="MonoDc")

DimPlot(Mono, reduction = "umap",label=T,repel = T)
ggsave('Monocyte_UMAP.pdf', dpi=900,width=7, height=5, units="in")



######## Fig 3E. DotPlot_(Monocyte_Cluster Markers)#########

Mono1<-Mono[,Idents(Mono) %in% c("cMono S100A8/9/12Hi","cMono MHCIIhi","NcMono","Type I IFN Mono","MonoDc")]
markers <- FindAllMarkers(Mono1, max.cells.per.ident = 200)

#write.csv(markers,"Mono_Marker_List.csv")


top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_markers<-unique(top10$gene)


  DotPlot(object = Mono1, assay="RNA", features  = top5_markers,cols="RdYlBu" ) + scale_y_discrete(limits = rev) +  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
     guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + guides(size = guide_legend(title = "Percent\nExpressed"),color = guide_colorbar(title = "Average\nExpression"))+ theme(axis.text.x  =element_text( hjust=1, size=10,angle=90))
ggsave('Monocyte_DotPlot.pdf', dpi=900,width=12, height=5, units="in")



######## Fig 3F. Volcano Plot (DEG, IPF vs HC in Monocytes)#########
# (code below is repeated for each ncMono subsets)

result<-read.delim("cMono_DEG_0.1.txt")
result$geneid<-rownames(result)


result$FDR <- ifelse(result$padj < 0.1 , "FDR < 0.1", "ns")
result$DEG <- ifelse(result$padj < 0.1 & abs(result$log2FoldChange) > 1.0 , "FDR < 0.1 and |Log2FC| 1.0", "ns")

result1<-result[!(is.na(result$DEG)),]


ggplot(result1, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point(aes(color = DEG))  +ggtitle("cMono: IPF vs HC")+
        scale_color_manual(values = c("red", "grey30")) + labs(x = "Log2 Fold Change") +
        theme(legend.position = "bottom",axis.text=element_text(size=13), plot.title = element_text(size = 15, face = "bold"))+
        geom_text_repel(
              data = subset(result, DEG =="FDR < 0.1 and |Log2FC| 1.0"),
              aes(label = geneid),
              size = 3.5,
              box.padding = unit(0.35, "lines"),
              point.padding = unit(0.3, "lines")
                )
ggsave('cMono_Volcano.pdf', dpi=900,width=7, height=5, units="in")


######## Fig 3G. Regulon Activity Heatmap of IFN TFs(Monocyte subset)#########

Blood<-readRDS("Blood_obj.rds")

Mono_b<-Blood[,Blood$CellType %in% ("cMono","ncMono")]


regulonAUC <- importAUCfromText("aucell_result.csv")
meta<-as.data.frame(Mono_b@meta.data)

meta$Group<-paste0(meta$CellType,"_",meta$Disease)


regulonActivity_byCellType <- sapply(split(rownames(meta), meta$Group),
                              function(cells) rowMeans(getAUC(regulonAUC)[,cells]))


head(regulonActivity_byCellType)

IFN_TF<-c(ELF1...,ATF3...,STAT2...,IRF8...,STAT1...,NFKB1...,IRF7...,USF1...,IRF2...,STAT3...,NR3C1...,ETV6...,KDM5A...,HIF1A...,RB1...,TCF4...,CHD1...,HCFC1...,IRF1...)
IFN_TF_f<-stringr::str_replace(IFN_TF, pattern="...$", "(+)")




TF_regulon<-regulonActivity_byCellType[rownames(regulonActivity_byCellType) %in% IFN_TF_f,]
TF_regulon_Scaled <- t(scale(t(TF_regulon), center = T, scale=T))



cairo_pdf("IFN_TFs_heatmap.pdf", height=8,width=5)
pheatmap::pheatmap(TF_regulon_Scaled,cluster_cols=FALSE,cluster_rows=TRUE, border_color = "black",fontsize_row=7.0)
dev.off()

######## Fig 3H. DotPlot of IFN TF regulon activity comparison #########

scenic <- read.csv("aucell1.csv")

rownames(scenic) <- scenic[, 1]
scenic <- t(as.matrix(scenic[, -1]))


Mono1[["SCENIC"]] <- CreateAssayObject(data = scenic[, Cells(Mono1)])

Mono1@active.assay <- "SCENIC"


Idents(Mono1)<-Mono1$CellType

for( cluster in unique(Mono1@active.ident)){
  print(cluster)

  mk1 <- FindMarkers(Mono1[, Mono1@active.ident == cluster], "IPF", "HC", group.by="Disease", logfc.threshold = 0, assay = "SCENIC", slot = "data",  test.use = "LR",
                                          max.cells.per.ident = 200)
  mk1$TF <- rownames(mk1)
  write.table(mk1, file=paste0(make.names(cluster), ".TF.network.cluster.diffsIPF_vs_HC.xls"), sep="\t", quote=FALSE, row.names = FALSE)

}




combn1<-read.csv("TF_FCandPadj_table_Mostafavi_Mono.csv")
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

ggsave('DotPlot_TFComparison_Mostafavi_Mono.pdf', dpi=900,width=4, height=6, units="in")
