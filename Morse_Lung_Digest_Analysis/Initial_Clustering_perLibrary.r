# Title     : SC14NOR
# Objective : SC14NOR sample initial QC and clustering
# Created by: JW
# Created on: 15/04/2024

library(Seurat)
library(ggplot2)
library(DropletUtils)
library(dplyr)
library(Matrix)
library(reshape2)
library(cluster)
library(fitdistrplus)

dat.dir <- "/ceph/project/holab/jwoo/IPF_HT/Integration/Morse/CR_processed_fromGEO/SC14NOR"
raw.mat <- Read10X(dat.dir)


cell.calls <- DropletUtils::emptyDrops(raw.mat)
head(cell.calls)
table(cell.calls$FDR < 0.05)

sample1 <- CreateSeuratObject(raw.mat, project = "SC14NOR", meta.data = as.data.frame(cbind(cell.calls)))
sample1 <- PercentageFeatureSet(sample1, pattern = "MT-", assay="RNA", col.name = "percent.mito")
sample1 <- PercentageFeatureSet(sample1, pattern = "RP", assay="RNA", col.name = "percent.rp")

#pre-filter
sample1 <- sample1[, sample1$FDR < 0.05 | sample1$nFeature_RNA > 100]
sample1
sample1@active.assay <- "RNA"
sample1 <- NormalizeData(sample1)
sample1 <- FindVariableFeatures(sample1)
sample1 <- ScaleData(sample1)
sample1 <- RunPCA(object = sample1, verbose = FALSE)
ElbowPlot(sample1, ndims = 50)
sample1 <- RunUMAP(object = sample1, dims = 1:20, verbose = FALSE)
sample1 <- FindNeighbors(object = sample1, dims = 1:20, verbose = FALSE)
sample1 <- FindClusters(object = sample1, verbose = FALSE, resolution = .7)
DimPlot(object = sample1, label = TRUE, reduction="umap")

qplot(sample1$nFeature_RNA, geom="density")
FeaturePlot(sample1, "nFeature_RNA")
FeaturePlot(sample1, "percent.mito")

#this will need further tweaking on cell type by cell type basis
sample1$Filt <- sample1$nFeature_RNA > 300 & sample1$percent.mito < 30
DimPlot(object = sample1, group.by="Filt", reduction="umap")
sample1 <- sample1[, sample1$Filt]
sample1

sample1 <- NormalizeData(sample1)
sample1 <- FindVariableFeatures(sample1)
sample1 <- ScaleData(sample1)
sample1 <- RunPCA(object = sample1, verbose = FALSE)
ElbowPlot(sample1, ndims = 50)
sample1 <- RunUMAP(object = sample1, dims = 1:20, verbose = FALSE)
sample1 <- FindNeighbors(object = sample1, dims = 1:20, verbose = FALSE)
sample1 <- FindClusters(object = sample1, verbose = FALSE, resolution = .7)
DimPlot(object = sample1, label = TRUE, reduction="umap")


sample1$Chemistry <- "V1"
sample1$Disease_Status <- "Normal_Control"
sample1$Disease <- "Normal"
sample1$Age <- 76
sample1$Sex <- "M"




saveRDS(sample1, "SC14NOR.RDS")
