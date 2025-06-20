#kalliso abundance files has been input as txi

library("tximport")
library("DESeq2")
library("AnnotationHub")
library("ensembldb")
library( "pheatmap" )
library( "RColorBrewer" )
library( "ggplot2" )
library( "biomaRt" )
library("dplyr")

colnames(txi$counts)

sampleTable <- data.frame(sample = colnames(df),
                          condition = c(rep("condition1", 3), rep("condition2", 3))) #same order when you input the files 

rownames(sampleTable) <- colnames(df)

dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ condition)

dds <- DESeqDataSetFromMatrix(countData = round(df), colData = sampleTable, design = ~ condition)

#Pre-filtering the dataset:
nrow(dds)
keep <- rowSums(counts(dds) >= 10) >=3
nrow(dds[keep,])
dds <- dds[keep,]
nrow(dds)

dds$condition <- relevel(dds$condition, ref = "condition2") #input the one be compared with 


#rlog transform for pca plot
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

sampleDists <- dist(t(assay(rld)))
sampleDists

#heatmap
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rownames(sampleDistMatrix), 
                                     sampleTable$condition, 
                                     sep = " - " )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
breaks <- seq(0, 80, length.out = 256)

p <- pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         breaks = breaks,
         show_rownames = TRUE,
         show_colnames = TRUE)

pdf(file.path(outdir, "heatmap_name.pdf"),
    height = 4,
    width = 5)
p
dev.off()


z <- plotPCA(rld, intgroup = c("condition"))

pdf(file.path(outdir, "pca_name.pdf"),
    height = 4,
    width = 5)
z + ggtitle(" PCA plot name")
dev.off()

#Results
dds <- DESeq(dds)

res <- results(dds, 
               contrast=c("condition","condition1","condition2"), #342 is the item that being compared
               alpha = 0.05,
               )

head(mcols(dds, use.names = TRUE))
summary(res, )

#MA plot
pdf(file.path(outdir, "MAplot_name.pdf"),
    height = 4,
    width = 5)
plotMA( res, ylim = c(-10, 10), alpha = 0.05 , main = "MA plot")
dev.off()



#Add external gene name

mart <- useEnsembl(biomart = 'ensembl', 
                   dataset = 'hsapiens_gene_ensembl')

genemap <- getBM( attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                  filters = c("ensembl_gene_id"),
                  values = rownames(res),
                  mart = mart )
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)
res_df <- res_df[,c(7,1,2,3,4,5,6)]

rownames(res_df) <- NULL

res_df <- merge(res_df, genemap, by = "ensembl_gene_id")
head(res_df)



#save deseq2 files
write.csv(res_df[order(res_df$padj),],
            file = file.path(outdir, "name_deseq2.csv"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

padj_0.05 <- subset(res_df, res_df$padj <= 0.05)

write.csv(padj_0.05[order(padj_0.05$log2FoldChange),],
            file = file.path(outdir, "name_deseq2_padj0.05.csv"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

log2foldchange_1_padj_0.05 <- subset(res_df, abs(res_df$log2FoldChange) >= 1 & res_df$padj <= 0.05)

write.csv(log2foldchange_1_padj_0.05[order(log2foldchange_1_padj_0.05$padj),],
            file = file.path(outdir, "name_deseq2_fc1_padj0.05.csv"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

log2foldchange_0.5_padj_0.05 <- subset(res_df, abs(res_df$log2FoldChange) >= 0.5 & res_df$padj <= 0.05)

write.csv(log2foldchange_0.5_padj_0.05[order(log2foldchange_0.5_padj_0.05$padj),],
            file = file.path(outdir, "name_deseq2_fc0.5_padj0.05.csv"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

