# Title     : Sample merge
# Objective : merge samples
# Created by: JW
# Created on: 05/02/2024

library(Seurat)
library(ggplot2)
library(DropletUtils)
library(dplyr)
library(Matrix)
library(reshape2)
library(cluster)
library(fitdistrplus)


samples <- lapply(list.files("/Path/to/individual_samples/",
                            pattern = ".RDS", full.names = T), readRDS)

merged <- merge(samples[[1]], samples[[2]])

##double check barcode overlap
for( i in 1:2){

 for( j in 1:2){
   if( i == j){
     next
   }else{

     print(table(Cells(samples[[i]]) %in% Cells(samples[[j]])))
   }
 }
}


setwd("/working_dir/")

###prelim clustering
merged@active.assay <- "RNA"
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(object = merged, verbose = FALSE)
ElbowPlot(merged, ndims = 50)
merged <- RunUMAP(object = merged, dims = 1:20, verbose = FALSE)
merged <- FindNeighbors(object = merged, dims = 1:20, verbose = FALSE)
merged <- FindClusters(object = merged, verbose = FALSE, resolution = .7)
DimPlot(object = merged, label = TRUE, reduction="umap")

DimPlot(object = merged, group.by="orig.ident", reduction="umap")


merged <- harmony::RunHarmony(merged, group.by.vars = "orig.ident")
merged <- RunUMAP(object = merged, dims = 1:20, verbose = FALSE, reduction = "harmony")

merged <- FindNeighbors(object = merged, dims = 1:20, verbose = FALSE,reduction = "harmony")
merged <- FindClusters(object = merged, verbose = FALSE, resolution = .7)
DimPlot(object = merged, label = TRUE, reduction="umap")


merged$Filt <- merged$nFeature_RNA > 500 & merged$percent.mito < 10 & merged$percent.rp > 5 |
 (merged$nFeature_RNA > 500 & merged$percent.mito < 30 & merged@active.ident %in% c("19"))

DimPlot(object = merged, group.by="Filt", reduction="umap")

merged <- merged[, merged$Filt]
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(object = merged, verbose = FALSE)
ElbowPlot(merged, ndims = 50)

merged <- harmony::RunHarmony(merged, group.by.vars = "orig.ident")
merged <- RunUMAP(object = merged, dims = 1:20, verbose = FALSE, reduction = "harmony")
merged <- FindNeighbors(object = merged, dims = 1:20, verbose = FALSE,reduction = "harmony")
merged <- FindClusters(object = merged, verbose = FALSE, resolution = .7)

DimPlot(object = merged, label = TRUE, reduction="umap")
DimPlot(object = merged, group.by="orig.ident", reduction="umap")
DimPlot(object = merged, group.by="Chemistry", reduction="umap")




markers <- FindAllMarkers(merged, max.cells.per.ident = 200)


merged <- RenameIdents(merged, "0"="FABP4 AM")
merged <- RenameIdents(merged, "2"="IGF1+ AM")
merged <- RenameIdents(merged, "3"="CD206hi FNhi AM")
merged <- RenameIdents(merged, "5"="monoSPP1+ AM")
merged <- RenameIdents(merged, "12"="CXCL10+ AMs")

merged <- RenameIdents(merged, "9"="Cycling Cells")
merged <- RenameIdents(merged, "10"="Cycling Cells")

merged <- RenameIdents(merged, "1"="CD4+ T-Cells")
merged <- RenameIdents(merged, "7"="CD4+ T-Cells")

merged <- RenameIdents(merged, "4"="CD8+ T-Cells")
merged <- RenameIdents(merged, "8"="NK Cells")

merged <- RenameIdents(merged, "6"="Macrophages/Monocytes")
merged <- RenameIdents(merged, "14"="Macrophages/Monocytes")
merged <- RenameIdents(merged, "20"="Macrophages/Monocytes")

merged <- RenameIdents(merged, "13"="B-Cells")

merged <- RenameIdents(merged, "15"="DCs")
merged <- RenameIdents(merged, "16"="DCs")
merged <- RenameIdents(merged, "18"="pDCs")

merged <- RenameIdents(merged, "19"="Epithelium")

merged <- RenameIdents(merged, "16"="DCs")

merged <- RenameIdents(merged, "17"="Mast")

saveRDS(merged, file="merged.RDS")

BAL<-merged[,merged$Tissue=="BAL"]
Blood<-merged[,merged$Tissue=="Blood"]


saveRDS(BAL, file="BAL_obj.RDS")
saveRDS(Blood, file="Blood_obj.RDS")
