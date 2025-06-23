library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(ReactomePA)
library(data.table)

line_of_interest = "3110_uqcrc2"
df <- read.csv("C:/Users/Jiawan/Dropbox/Jiawan/Analysis/Deseq2/comparison/CRISPRa_fam3_new_result/3110_uqcrc2_versus_3110_ev/npc_3110_uqcrc2_versus_3110_ev_deseq2_fc0.5_padj0.05.csv")
#C:/Users/Jiawan/Dropbox/Jiawan/Analysis/Deseq2/comparison/results/CRISPRa/reverse_DEGs

outdir <- "C:/Users/Jiawan/Dropbox/Jiawan/Analysis/CRISPRa_GO_term/gene"

df[which(df$log2FoldChange >= 0.5),'sig'] <- 'up'
df[which(df$log2FoldChange <= -0.5),'sig'] <- 'down'

df_up <- subset(df, sig == 'up')
df_down <- subset(df, sig == 'down')

diff_up <- df_up 

gene_up.df <- bitr(diff_up$ensembl_gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) #change$SYMBOL based on column name
gene_up <- gene_up.df$ENTREZID

ego_ALL <- enrichGO(gene = gene_up,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   minGSSize = 5,
                   maxGSSize = 5000,
                   pvalueCutoff = 0.05, #could be 0.01
                   qvalueCutoff = 0.05,
                   readable = TRUE)


goplot <- dotplot(ego_ALL, showCategory = 20, font.size = 8) +
          ggtitle(paste("Dotplot for GO Enrichment -", line_of_interest, "_upregulated Genes"))

output_file_base <- paste0("npc_", line_of_interest, "_upregulated_genes")



write.csv(as.data.frame(ego_ALL),
          file = file.path(outdir, paste0("ego_ALL_", output_file_base, ".csv")),
          quote = FALSE,
          row.names = TRUE,
          col.names = TRUE)


diff_down <- df_down

gene_down.df <- bitr(diff_down$ensembl_gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) #change$SYMBOL based on column name
gene_down <- gene_down.df$ENTREZID

ego_ALL <- enrichGO(gene = gene_down,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   minGSSize = 5,
                   maxGSSize = 5000,
                   pvalueCutoff = 0.05, #could be 0.01
                   qvalueCutoff = 0.05,
                   readable = TRUE)


goplot <- dotplot(ego_ALL, showCategory = 20, font.size = 8) +
          ggtitle(paste("Dotplot for GO Enrichment -", line_of_interest, "_downregulated Genes"))

output_file_base <- paste0("npc_", line_of_interest, "_downregulated_genes")



write.csv(as.data.frame(ego_ALL),
          file = file.path(outdir, paste0("ego_ALL_", output_file_base, ".csv")),
          quote = FALSE,
          row.names = TRUE,
          col.names = TRUE)


#KEGG
kk <- enrichKEGG(gene         = gene_up,
                 organism     = 'hsa',
                 qvalueCutoff = 0.05)
head(kk)

goplot <- dotplot(kk, showCategory=20, font.size=8) + ggtitle(paste("Dotplot for KEGG Enrichment -", line_of_interest, "_upregualted_Genes"))

output_file_base <- paste0("npc_", line_of_interest, "_upregulated_genes")



write.csv(as.data.frame(kk),
          file = file.path(outdir,  paste0("KEGG_", output_file_base, ".csv")),
          quote = FALSE,
          row.names = TRUE,
          col.names = TRUE)


kk <- enrichKEGG(gene         = gene_down,
                 organism     = 'hsa',
                 qvalueCutoff = 0.05)
head(kk)

goplot <- dotplot(kk, showCategory=20, font.size=8) + ggtitle(paste("Dotplot for KEGG Enrichment -", line_of_interest, "_downregualted_Genes"))

output_file_base <- paste0("npc_", line_of_interest, "_downregulated_genes")



write.csv(as.data.frame(kk),
          file = file.path(outdir,  paste0("KEGG_", output_file_base, ".csv")),
          quote = FALSE,
          row.names = TRUE,
          col.names = TRUE)

#reactome
pathway.reac <- enrichPathway(gene_up)

head(pathway.reac)

goplot <- dotplot(pathway.reac, showCategory=20, font.size=8) + ggtitle(paste("Dotplot for Reactome Enrichment -", line_of_interest, "_upregulated_Genes"))

output_file_base <- paste0("npc_", line_of_interest, "_upregulated_genes")



write.csv(as.data.frame(pathway.reac),
          file = file.path(outdir, paste0("Reactome_", output_file_base, ".csv")),
          quote = FALSE,
          row.names = TRUE,
          col.names = TRUE)


pathway.reac <- enrichPathway(gene_down)

head(pathway.reac)

goplot <- dotplot(pathway.reac, showCategory=20, font.size=8) + ggtitle(paste("Dotplot for Reactome Enrichment -", line_of_interest, "_downregulated_Genes"))

output_file_base <- paste0("npc_", line_of_interest, "_downregulated_genes")



write.csv(as.data.frame(pathway.reac),
          file = file.path(outdir, paste0("Reactome_", output_file_base, ".csv")),
          quote = FALSE,
          row.names = TRUE,
          col.names = TRUE)
