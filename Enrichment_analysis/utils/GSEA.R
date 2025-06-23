require(tidyverse)
require(ggplot2)
require(clusterProfiler)
require(org.Hs.eg.db)
require(enrichplot)

#input deseq2 file:the one for volacano plot

DEG_3<- fam3_im_3110_uqcrc2_331_vs_321_ev_deseq2_padj0.05
head(DEG_3)

genelist<-DEG_3$external_gene_name
genelist<- genelist %>%  bitr(fromType = "SYMBOL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)

colnames(genelist)[1]<-"genename"
colnames(DEG_3)[8]<-"genename"

DEG_3<-genelist %>%
  inner_join(DEG_3,by='genename') %>% 
  ## select配合everything排序把，改变变量顺序
  dplyr::select(ENTREZID,log2FoldChange,everything())
head(DEG_3)

geneList<-DEG_3$log2FoldChange
names(geneList)<-as.character(DEG_3$ENTREZID)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
head(geneList)

require(enrichplot)
library(msigdbr)


m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene)

head(m_t2g)

em <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2g,pvalueCutoff = 0.1,verbose = FALSE)
em[1:5,1:5]

write.csv(em,file = "npc_3110_uqcrc2_versus_ev_controls_GSEA_reactome.csv",row.names = T)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>% 
  dplyr::select(gs_name, entrez_gene)

em_1 <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2g,pvalueCutoff = 0.1,verbose = FALSE)
write.csv(em_1,file = "npc_3110_uqcrc2_versus_ev_controls_GSEA_wikipathway.csv",row.names = T)

m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

em <- GSEA(geneList,nPerm = 10000 ,TERM2GENE = m_t2g,pvalueCutoff = 0.1,verbose = FALSE)
write.csv(em,file = "npc_3110_uqcrc2_versus_ev_controls_GSEA.csv",row.names = T)
