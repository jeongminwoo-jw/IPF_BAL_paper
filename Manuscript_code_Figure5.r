# Title     : LRT analysis_and_Fig5E
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


subs<-morse_sub[,morse_sub$CellType_v2 == Celltype[i]]

mtx_aggr<-AggregateExpression(subs,assays = "RNA",slot="count",group.by = "Donor")
cluster_counts<-as.data.frame(mtx_aggr$RNA)


################# Metadata generation ################

 metadata<-as.data.frame(subs@meta.data)

 meta_1<-metadata[,c("Donor","Disease_1")]

 meta_f<-unique(meta_1)

 rownames(meta_f)<-meta_f$Donor

meta_order<-meta_f[order(meta_f$Disease_1),]

 cluster_counts_f<-cluster_counts[,rownames(meta_order)]

meta_final<-meta_order
 ########### DESeq2 Run #########################

 meta_final$Sample_Name<-as.factor(meta_final$Donor)
 meta_final$Status<-as.factor(meta_final$Disease_1)

 dds <- DESeqDataSetFromMatrix(cluster_counts_f,
                               colData = meta_final,
                               design = ~ Status)


 #Prefilter
 dds <- dds[rowSums(counts(dds))> 3, ]
 rld <- rlog(dds, blind=TRUE)

LRT_DEG<-read.csv(paste0(Celltype[i],"_Cluster_annotation.csv"))

 DEG_mat<-as.data.frame(assay(rld)[rownames(assay(rld)) %in% LRT_DEG$X,])

DEG_mat$IPF_L_median = rowMedians(as.matrix(DEG_mat[,1:5]))
DEG_mat$IPF_U_median = rowMedians(as.matrix(DEG_mat[,6:8]))
DEG_mat$HC_median = rowMedians(as.matrix(DEG_mat[,9:12]))

DEG_mat$Cluster<-LRT_DEG$Cluster

write.csv(DEG_mat,paste0(Celltype[i],"_Normalised_Expression_mat.csv"))


DEG_mat_1<-DEG_mat[,13:15]
DEG_mat_1$gene<-rownames(DEG_mat_1)
data1<-melt(DEG_mat_1)


data1$Cluster<-ifelse((data1$gene) %in% LRT_DEG[LRT_DEG$Cluster =="1",]$X, "1","2")



data1$Cluster<-data1$gene

data1$Cluster<-ifelse((data1$gene) %in% LRT_DEG[LRT_DEG$Cluster =="1",]$X, "1",data1$gene)
data1$Cluster<-ifelse((data1$gene) %in% LRT_DEG[LRT_DEG$Cluster =="2",]$X, "2",data1$Cluster)
data1$Cluster<-ifelse((data1$gene) %in% LRT_DEG[LRT_DEG$Cluster =="3",]$X, "3",data1$Cluster)
data1$Cluster<-ifelse((data1$gene) %in% LRT_DEG[LRT_DEG$Cluster =="4",]$X, "4",data1$Cluster)



########## Fig 5E BarPlot ######################

p<-ggplot(data1, aes(variable,value, fill=Cluster)) +
     geom_boxplot()+

     # geom_line() joins the paired datapoints
     # color and size parameters are used to customize line
     geom_line(aes(group = gene), size=0.3, color='gray60', alpha=0.4)+

     # geom_point() is used to make points at data values
     # fill and size parameters are used to customize point
     geom_point(aes(fill=Cluster,group=gene),size=3,shape=21)


cairo_pdf(paste0(Celltype[i],"_boxPlot.pdf"),width=7, height=5)
print(p)
dev.off()


}


########## Fig 5F. Pathway analysis ######################

# (Repeated for each clusters for Celltypes)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggnewscale)
library(fgsea)
library(dplyr)
library(msigdbr)

Gene_anno<-read.csv("Cluster_annotation.csv")

Gene_GO<-Gene_anno[Gene_anno$Cluster =="1",]
gene_up<-as.vector(Gene_GO$X)


########## GO  #########


ego_up<-enrichGO(gene_up,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont="all",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)

head(summary(ego_up))

ego_up

cluster_summary<-data.frame(ego_up)

write.table(cluster_summary,"Cluster1_GO.txt",sep="\t")

dotplot(ego_up, split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale="free") + scale_fill_viridis(direction = -1)



ego_BP<-enrichGO(gene_up, OrgDb= org.Mm.eg.db,keyType = "SYMBOL",ont="BP",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
dotplot(ego_BP,showCategory=20,label_format=100)+ scale_fill_viridis(direction = -1)

ego2_BP <- simplify(ego_BP)
dotplot(ego2_BP,showCategory=20,label_format=100)+ scale_fill_viridis(direction = -1)

cnetplot(ego2_BP, colorEdge = TRUE)


########## Reactome #########

gmtfile<-"h.all.v6.1.symbols.gmt"
Hallmark<-read.gmt(gmtfile)
enrich_Hallmark<-enricher(gene,TERM2GENE = Hallmark,minGSSize = 1,maxGSSize = 1000)

cluster_summary<-data.frame(enrich_Hallmark)
write.table(cluster_summary,"Cluster1_enrichment.txt",sep="\t")
dotplot(enrich_Hallmark,showCategory=32)
