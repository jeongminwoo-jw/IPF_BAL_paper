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

merged <- merge(samples[[1]], samples[2:16])

##double check barcode overlap
for( i in 1:16){

 for( j in 1:16){
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
DimPlot(object = merged, group.by="Chemistry", reduction="umap")


merged <- harmony::RunHarmony(merged, group.by.vars = "orig.ident")
merged <- RunUMAP(object = merged, dims = 1:20, verbose = FALSE, reduction = "harmony")

merged <- FindNeighbors(object = merged, dims = 1:20, verbose = FALSE,reduction = "harmony")
merged <- FindClusters(object = merged, verbose = FALSE, resolution = .7)
DimPlot(object = merged, label = TRUE, reduction="umap")


merged$Filt <- merged$nFeature_RNA > 500 & merged$percent.mito < 10 & merged$percent.rp > 5 |
 (merged$nFeature_RNA > 500 & merged$percent.mito < 30 & merged@active.ident %in% c("0", "1","3"))

DimPlot(object = merged, group.by="Filt", reduction="umap")

mk <- FindMarkers(merged, "19", max.cells.per.ident = 200)

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

saveRDS(merged, file="merged.RDS")
