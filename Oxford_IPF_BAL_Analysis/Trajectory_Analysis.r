
library(Seurat)
library(SingleCellExperiment)
library(scales)
library(SummarizedExperiment)
library(slingshot)
library(viridis)
library(tradeSeq)
library(UpSetR)

################################################################################
############### Trajectory analysis using slingshot ############################
################################################################################


setwd("/ceph/project/holab/jwoo/IPF_HT/Trajectory/slingshot/Monocyte_Root_woCyclingCells")

######## Load seurat object of Mono_Mac subset

MacroMono_obj<-readRDS("../Monocle3/MacroMono_obj.rds")


######## Assign monocyte as a Root ##########################
MacroMono_obj$Celltype_1<-MacroMono_obj$Celltype
MacroMono_obj$Celltype_1<-gsub("CD14[+] Monocytes","Monocytes",MacroMono_obj$Celltype_1)
MacroMono_obj$Celltype_1<-gsub("CD16[+] Monocytes","Monocytes",MacroMono_obj$Celltype_1)



DefaultAssay(object = MacroMono_obj) <- "RNA"


sds <- slingshot(Embeddings(MacroMono_obj, "UMAP"), clusterLabels = MacroMono_obj$Celltype_1,
                start.clus = c("Monocytes"))

SCE <- as.SingleCellExperiment(MacroMono_obj)




#' Assign a color to each cell based on some value
#'
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.

cell_pal <- function(cell_vars, pal_fun,...) {
            if (is.numeric(cell_vars)) {
            pal <- pal_fun(100, ...)
            return(pal[cut(cell_vars, breaks = 100)])
            } else {
            categories <- sort(unique(cell_vars))
            pal <- setNames(pal_fun(length(categories), ...), categories)
            return(pal[cell_vars])
            }
        }

cell_colors_clust <- cell_pal(MacroMono_obj$Celltype, hue_pal())

plot(reducedDim(SCE,"UMAP"), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(SlingshotDataSet(sds), lwd = 2, type = 'lineages', col = 'black')

nc <- 2
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)

pal <- viridis(100, end = 0.95)


par(mar=c(1,1,1,1))
par(mfrow = c(nr, nc))
for (i in nms) {
     colors <- pal[cut(pt[,i], breaks = 100)]
     plot(reducedDim(SCE,"UMAP"), col = colors, pch = 16, cex = 0.5, main = i)
     lines(SlingshotDataSet(sds), lwd = 2, col = 'black', type = 'lineages')
 }


 write.csv(pt,"SlingShot_Pseudotime.csv")
 saveRDS(sds,"Merged_SlingShot.RDS")
# Plot reformat
colors <- pal[cut(pt[,1], breaks = 100)]
plot(reducedDim(SCE,"UMAP"), col = colors, pch = 16, cex = 0.5, main = "Lineage1")
lines(SlingshotDataSet(sds), lwd = 2, col = 'black', type = 'lineages')


 # legend
lgd <- matrix(hcl.colors(50), nrow=1)
rasterImage(lgd, -5,-7,-2,-6.7)
text(c(-5,-2), c(-6.8,-6.8), pos = 1, cex = .7,
       labels = format(range(slingPseudotime(sds)[,1], na.rm = TRUE), digits = 3))
text(-3.5, -6.9, pos = 3, cex = .7, labels = 'Pseudotime')


colors <- pal[cut(pt[,2], breaks = 100)]
plot(reducedDim(SCE,"UMAP"), col = colors, pch = 16, cex = 0.5, main = "Lineage2")
lines(SlingshotDataSet(sds), lwd = 2, col = 'black', type = 'lineages')

 # legend
 lgd <- matrix(hcl.colors(50), nrow=1)
 rasterImage(lgd, -5,-7,-2,-6.7)
 text(c(-5,-2), c(-6.8,-6.8), pos = 1, cex = .7,
       labels = format(range(slingPseudotime(sds)[,2], na.rm = TRUE), digits = 3))
 text(-3.5, -6.9, pos = 3, cex = .7, labels = 'Pseudotime')



 #### Density Plot   ###########

pt<-as.data.frame(pt)
pt$barcode<-rownames(pt)

 Meta<-as.data.frame(MacroMono_obj@meta.data)
 Meta$barcode<-rownames(Meta)
 Meta_1<-Meta %>% select(Celltype,barcode)

 slingshot_df<-pt %>% left_join(Meta_1,by="barcode")
 rownames(slingshot_df)<-slingshot_df$barcode


 ggplot(slingshot_df, aes(Lineage1, fill = Celltype, colour = Celltype) )+
    geom_density(alpha = 0.1) + theme_bw()





################################################################################
############### Association test using tradeSeq ############################
################################################################################


#### fitting the NB-GAM model  (Run as background job )####

set.seed(3)
sce <- fitGAM(sds1, conditions = factor(sds1@colData$Disease),
               nknots = 5)


saveRDS(sce,"fitGAM_result.rds")


####### Association test after NB-GAM fitting   #########


sce<-readRDS("fitGAM_result.rds")

rowData(sce)$assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce)$assocRes


head(assocRes)

# Lineage 1
IPF_Genes <-  rownames(assocRes)[
     which(p.adjust(assocRes$pvalue_lineage1_conditionIPF, "fdr") <= 0.05)]
HC_Genes <-  rownames(assocRes)[
     which(p.adjust(assocRes$pvalue_lineage1_conditionHC, "fdr") <= 0.05)]

UpSetR::upset(fromList(list(IPF = IPF_Genes, HC = HC_Genes)))


# Lineage 2
IPF_Genes <-  rownames(assocRes)[
     which(p.adjust(assocRes$pvalue_lineage2_conditionIPF, "fdr") <= 0.05)]
HC_Genes <-  rownames(assocRes)[
     which(p.adjust(assocRes$pvalue_lineage2_conditionHC, "fdr") <= 0.05)]

UpSetR::upset(fromList(list(IPF = IPF_Genes, HC = HC_Genes)))

write.csv(assocRes,"Associate_Test.csv")




################## Condition test ##################

condRes <- conditionTest(sce, lineages=T, l2fc = log2(2))

head(condRes)


condRes$padj_global <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)

sum(condRes$padj <= 0.1, na.rm = TRUE)


conditionGenes <- rownames(condRes)[condRes$padj <= 0.1]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]


###### Lineage1 DEG


condRes$padj_lineage1 <- p.adjust(condRes$pvalue_lineage1, "fdr")
mean(condRes$padj_lineage1 <= 0.05, na.rm = TRUE)

sum(condRes$padj_lineage1 <= 0.1, na.rm = TRUE)


conditionGenes_L1 <- rownames(condRes)[condRes$padj_lineage1 <= 0.1]
conditionGenes_L1 <- conditionGenes_L1[!is.na(conditionGenes_L1)]



###### Lineage2 DEG


condRes$padj_lineage2 <- p.adjust(condRes$pvalue_lineage2, "fdr")
mean(condRes$padj_lineage2 <= 0.05, na.rm = TRUE)

sum(condRes$padj_lineage2 <= 0.1, na.rm = TRUE)


conditionGenes_L2 <- rownames(condRes)[condRes$padj_lineage2 <= 0.1]
conditionGenes_L2 <- conditionGenes_L2[!is.na(conditionGenes_L2)]

write.csv(condRes,"Condition_Test.csv")



##########  Heatmaps of genes DE between conditions



### based on mean smoother
yhatSmooth <- predictSmooth(sce, gene = conditionGenes_L1, nPoints = 50, tidy = FALSE)


######### Fetch Lineage 1 smoother
yhatSmooth_1<-yhatSmooth [, c(1:50,101:150)]

yhatSmoothScaled <- t(scale(t(yhatSmooth_1)))
heatSmooth_HC_L1 <- pheatmap(yhatSmoothScaled[, 1:50],
cluster_cols = FALSE,
show_rownames = FALSE, show_colnames = FALSE, main = "HC Lineage1", legend = FALSE,
silent = TRUE
)

matchingHeatmap_IPF_L1 <- pheatmap(yhatSmoothScaled[heatSmooth_HC_L1$tree_row$order, 51:100],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = FALSE, show_colnames = FALSE, main = "IPF Lineage1",
  legend = FALSE, silent = TRUE
)

grid.arrange(heatSmooth_HC_L1[[4]], matchingHeatmap_IPF_L1[[4]], ncol = 2)




######### Fetch Lineage 2 smoother


yhatSmooth <- predictSmooth(sce, gene = conditionGenes_L2, nPoints = 50, tidy = FALSE)

yhatSmooth_1<-yhatSmooth[, c(51:100,151:200)]

yhatSmoothScaled <- t(scale(t(yhatSmooth_1)))
heatSmooth_HC_L2 <- pheatmap(yhatSmoothScaled[, 1:50],
cluster_cols = FALSE,
show_rownames = T, show_colnames = FALSE, main = "HC Lineage2", legend = FALSE,
silent = TRUE
)

matchingHeatmap_IPF_L2 <- pheatmap(yhatSmoothScaled[heatSmooth_HC_L2$tree_row$order, 51:100],
  cluster_cols = FALSE, cluster_rows = FALSE,
  show_rownames = T, show_colnames = FALSE, main = "IPF Lineage2",
  legend = FALSE, silent = TRUE
)

grid.arrange(heatSmooth_HC_L2[[4]], matchingHeatmap_IPF_L2[[4]], ncol = 2)



#######  Clustering gene sets with co-expression within H-clustering



row_anno <-as.data.frame(cutree(heatSmooth_HC_L1$tree_row, 7))
colnames(row_anno) <- "Cluster"
row_anno$Cluster <- as.factor(row_anno$Cluster)



yhatSmooth <- predictSmooth(sce, gene = conditionGenes_L1, nPoints = 50, tidy = FALSE)
yhatSmooth_1<-yhatSmooth [, c(1:50,101:150)]
yhatSmoothScaled <- t(scale(t(yhatSmooth_1)))


heatSmooth_HC_L1_1 <- pheatmap(yhatSmoothScaled[, 1:50],
cluster_cols = FALSE,
show_rownames = FALSE, show_colnames = FALSE, main = "A Lineage1", legend = FALSE,
silent = TRUE,annotation_row = row_anno
)

############ All Samples in one Heatmap ############


heatSmooth_HC_L1_2 <- pheatmap(yhatSmoothScaled[, 1:100],
cluster_cols = FALSE,
show_rownames = FALSE, show_colnames = FALSE, main = "Lineage1", legend = FALSE,
silent = TRUE
)

row_anno <-as.data.frame(cutree(heatSmooth_HC_L1_2$tree_row, 8))
colnames(row_anno) <- "Cluster"
row_anno$Cluster <- as.factor(row_anno$Cluster)
dev.off()


heatSmooth_HC_L1_2 <- pheatmap(yhatSmoothScaled[, 1:100],
cluster_cols = FALSE,
show_rownames = FALSE, show_colnames = FALSE, main = "Lineage1", legend = FALSE,
silent = TRUE,annotation_row = row_anno
)

write.csv(row_anno,"Cluster_Annotation_Lineage1.csv")


############ All Samples in one Heatmap (Lineage2) ############

yhatSmooth <- predictSmooth(sce, gene = conditionGenes_L2, nPoints = 50, tidy = FALSE)
yhatSmooth_1<-yhatSmooth[, c(51:100,151:200)]
yhatSmoothScaled <- t(scale(t(yhatSmooth_1)))


heatSmooth_HC_L2_2 <- pheatmap(yhatSmoothScaled[, 1:100],
cluster_cols = FALSE,
show_rownames = FALSE, show_colnames = FALSE, main = "Lineage2", legend = FALSE,
silent = TRUE
)


row_anno <-as.data.frame(cutree(heatSmooth_HC_L2_2$tree_row, 8))
colnames(row_anno) <- "Cluster"
row_anno$Cluster <- as.factor(row_anno$Cluster)
dev.off()

heatSmooth_HC_L2_2 <- pheatmap(yhatSmoothScaled[, 1:100],
cluster_cols = FALSE,
show_rownames = FALSE, show_colnames = FALSE, main = "Lineage2", legend = FALSE,
silent = TRUE,annotation_row = row_anno
)

write.csv(row_anno,"Cluster_Annotation_Lineage2.csv")



####### different EndTest

endRes <- diffEndTest(sce,l2fc = log2(2))
endRes$padj_global <- p.adjust(endRes$pvalue, "fdr")
head(endRes)

table(endRes$padj_global <= 0.1)

write.csv(endRes,"endTest_Lineages.csv")

endRes_gene <- rownames(endRes)[endRes$padj_global <= 0.1]
endRes_gene  <- endRes_gene [!is.na(endRes_gene )]


# Predict smooth

endSmooth <- predictSmooth(sce,gene = endRes_gene , nPoints = 50, tidy = FALSE)
endSmoothScaled <- t(scale(t(endSmooth)))



endTest_L1 <- pheatmap(endSmoothScaled,
cluster_cols = FALSE,
show_rownames = T, show_colnames = FALSE, main = "Differentiated Cell type Markers",
silent = TRUE
)



#####filtering out post bifurcation ########



sce<-readRDS("fitGAM_result.rds")
meta_slingshot<-read.csv("../Metadata_with_ModuleScore&Pseudotime.csv")

meta_L1<-meta_slingshot[meta_slingshot$Lineage1 >12.5,]
meta_L2<-meta_slingshot[meta_slingshot$Lineage2 >12.5,]

cells_keep<-unique(c(meta_L1$X,meta_L2$X))

sce1<-sce[,rownames(sce@colData) %in% cells_keep]
