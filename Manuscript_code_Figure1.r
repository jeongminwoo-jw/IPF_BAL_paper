require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(Hmisc)
library(RColorBrewer)




######## Fig 1B. UMAP #########

BAL<-readRDS("BAL_obj.rds")

print(head(BAL@meta.data))
DefaultAssay(object = BAL) <- "RNA"
Idents(BAL)<-BAL@meta.data$CellType_v2
DimPlot(BAL, reduction = "umap", group.by="Celltype_v2",cols=c("#d8eab2","grey","grey52","#D62728","#E377C2","#9467BD","#FF9896","#C49C94","#2CA02C" ,"#BCBD22","#8C564B", "#F7B6D2", "#C5B0D5", "#FFBB78", "#AEC7E8" , "#98DF8A", "#1F77B4"),label=T,repel = T)
ggsave('IPF_BAL_UMAP.pdf', dpi=900,width=7, height=5, units="in")


######## Fig 1D. DotPlot #########


BAL1<-BAL[,Idents(BAL) %in% c("FABP4 AM","IGF1+ AM","CD206hi FNhi AM","monoSPP1+ AM","CXCL10+ AMs")]
markers <- FindAllMarkers(BAL1, max.cells.per.ident = 200)


top10<-markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_markers<-unique(top10$gene)

Idents(BAL1)<-factor(BAL1$CellType_v2,levels = c("FABP4 AM","IGF1+ AM","CD206hi FNhi AM","CXCL10+ AMs","monoSPP1+ AM"))
DotPlot(object = BAL1, assay="RNA", features  = top5_markers,cols="RdYlBu" ) + scale_y_discrete(limits = rev) +  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + guides(size = guide_legend(title = "Percent\nExpressed"),color = guide_colorbar(title = "Average\nExpression"))+ theme(axis.text.x  =element_text( hjust=1, size=10,angle=90))

ggsave('Mac_Marker_DotPlot.pdf', dpi=900,width=10, height=4.1, units="in")



######## Fig 1F. TF motif DotPlot  #########

Idents(BAL)<-BAL@meta.data$CellType_v2

BAL@active.assay <- "SCENIC"


TF_markers<-read.delim("TF.network.Celltype.markers.xls")

TF_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

 DotPlot(merged, features=unique(top10$gene), assay = "SCENIC") + scale_colour_viridis_c(option = "D") + theme(axis.text.x = element_text(angle = 45, hjust=1))

 ggsave('TF_Celltype_Marker_DotPlot.pdf', dpi=900,width=11, height=4.5, units="in")



######## Fig 1G. DimPlot (Morse Lung data)  #########


Morse<-readRDS("Morse_merged.RDS")

DimPlot(Morse,group.by="Celltype",cols= c("#A6CDE2","#1E78B4","#74C476","#34A047","#F59899","#E11E26",
               "#FCBF6E","#F47E1F","#CAB2D6","#6A3E98","#FAF39B","#B15928"))

ggsave('MorseLung_UMAP.pdf', dpi=900,width=8.5, height=5, units="in")


######## Fig 1I. Celltype Marker DotPlot (Morse) #########



Idents(Morse)<-sort(Morse$CellType_v3)


markers <- FindAllMarkers(Morse, max.cells.per.ident = 200)

#write.csv(markers,"5_6_Marker_List.csv")


top5<-markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_markers<-unique(top5$gene)


DotPlot(object = Morse, assay="RNA", features  = top5_markers,cols="RdYlBu" ) + scale_y_discrete(limits = rev) +  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
     guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + guides(size = guide_legend(title = "Percent\nExpressed"),color = guide_colorbar(title = "Average\nExpression"))+ theme(axis.text.x  =element_text( hjust=1, size=10,angle=90))

ggsave('Morse_CelltypeMarker_DotPlot.pdf', dpi=900,width=8.5, height=5, units="in")



######## Fig 1H. DimPlot (Morse Myeloid Subsets)  #########
Morse_sub<-readRDS("Morse_Myeloid_MAY2024_f.rds")
 head(Morse_sub@meta.data)

DimPlot(Morse_sub, reduction = "umap", group.by = "CellType")
ggsave('MorseMyeloid_UMAP.pdf', dpi=900,width=y, height=5, units="in")
