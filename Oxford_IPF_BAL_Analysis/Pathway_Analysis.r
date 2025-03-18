library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggnewscale)
library(fgsea)
library(dplyr)
library(msigdbr)



##########################################################
############### Pathway analysis ##########################
###########################################################



DEG_results <-read.delim("AM1_DEG.txt",header=T)
DEGs <- DEG_results[!is.na(DEG_results$padj) & DEG_results$padj < 0.05,]
DEGS_up<-DEGs[abs(DEGs$log2FoldChange) > 1.0 ,]
DEGS_down<-DEGs[abs(DEGs$log2FoldChange) < -1.0 ,]

gene_up<-as.vector(DEGs_up$gene_id)


########## GO  #########


ego_up<-enrichGO(gene_up,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont="all",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)

head(summary(ego_up))

ego_up

cluster_summary<-data.frame(ego_up)

write.table(cluster_summary,"AM1_GO.txt",sep="\t")

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
write.table(cluster_summary,"AM1_enrichment.txt",sep="\t")
dotplot(enrich_Hallmark,showCategory=32)


##########################################################
############### GSEA analysis ############################
###########################################################


####Get Ranked Input from DEG Results

ranks <- DEG_results$log2FoldChange
names(ranks)<-rownames(DEG_results)
head(ranks)




#########Prepare GeneSet

## Hallmark
m_df_h<- msigdbr(species = "Homo sapiens", category = "H")
head(m_df_h)
fgsea_sets<- m_df_h %>% split(x = .$gene_symbol, f = .$gs_name)
head(fgsea_sets)

## GO
m_df_GO<- msigdbr(species = "Homo sapiens", category = "C5",subcategory ="GO:BP")
head(m_df_GO)
fgsea_sets_GO<- m_df_GO %>% split(x = .$gene_symbol, f = .$gs_name)
head(fgsea_sets_GO)

## Reactome
m_df_Reactome<- msigdbr(species = "Homo sapiens", category = "C2",subcategory ="CP:REACTOME")
head(m_df_Reactome)
fgsea_sets_Reactome<- m_df_Reactome %>% split(x = .$gene_symbol, f = .$gs_name)
head(fgsea_sets_Reactome)


#### HallMark_GeneSets

fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks,
                  nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>%
  head()


#plot a barplot for with the normalized Enrichment score


ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +   labs(x="Pathway", y="Normalized Enrichment Score", title="Hallmark pathways NES from GSEA")+
  theme_bw()


data.table::fwrite(fgseaRes, file = paste0("AM1", '_gsea_Hallmark.tsv'), sep = "\t", sep2 = c("", " ", ""))

#### GO terms

fgseaRes_GO <- fgsea(pathways = fgsea_sets_GO,
                  stats = ranks,
                  nperm=1000)


fgseaResTidy_GO <- fgseaRes_GO %>%
                    as_tibble() %>%
                    arrange(desc(NES))

fgseaResTidy_GO %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>%
  head()


# only plot the top 25 pathways

fgsea_go_sig<-fgseaResTidy_GO %>% filter(padj < 0.05)
fgsea_up<-fgsea_go_sig%>% head(n= 25)
fgsea_down<-fgsea_go_sig%>% tail(n= 25)

fgsea_up$Positive<-"TRUE"
fgsea_down$Positive<-"FALSE"

fgsea_go_final<-rbind(fgsea_up,fgsea_down)

fgsea_go_final$pathway<-gsub("GOBP_","",fgsea_go_final$pathway)

ggplot(fgsea_go_final, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=Positive)) +
    coord_flip() +    labs(x="Pathway", y="Normalized Enrichment Score",title="GO BiologicalProcess  NES from GSEA")+
    theme_bw()


data.table::fwrite(fgseaRes_GO , file = paste0("Macs", '_gsea_GO.tsv'), sep = "\t", sep2 = c("", " ", ""))



##REACTOME

fgseaRes_Reactome <- fgsea(pathways = fgsea_sets_Reactome,
                    stats = ranks,
                    nperm=1000)


fgseaResTidy_Reactome <- fgseaRes_Reactome %>%
                   as_tibble() %>%
                  arrange(desc(NES))

fgseaResTidy_Reactome%>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    arrange(padj) %>%
    head()


# only plot the top 25 pathways

fgsea_c2_sig<-fgseaResTidy_Reactome %>% filter(padj < 0.06)
fgsea_up<-fgsea_c2_sig%>% head(n= 25)
fgsea_down<-fgsea_c2_sig%>% tail(n= 25)

fgsea_up$Positive<-"TRUE"
fgsea_down$Positive<-"FALSE"

fgsea_c2_final<-rbind(fgsea_up,fgsea_down)

fgsea_c2_final$pathway<-gsub("REACTOME_","",fgsea_c2_final$pathway)

ggplot(fgsea_c2_final, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Positive)) +
  coord_flip() + labs(x="Pathway", y="NES", title="Reactome  NES from GSEA")+
  theme_bw()
