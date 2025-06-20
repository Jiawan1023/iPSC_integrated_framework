library(DESeq2)
library(ggplot2)

cnts <- iPSCsignal #input peak counts

rownames(cnts) = paste(cnts$chr, cnts$start, cnts$end, sep=":")

cnt1 = cnts[,-(1:4)]
# alternative way
rownames(cnts) <- paste(cnts$chr, cnts$start, cnts$end, sep = ":")
rownames(cnts) <- gsub(":", "-", rownames(cnts), fixed = TRUE)
rownames(cnts) <- sub("-", ":", rownames(cnts))  # Change first "-" to ":"

# change select column for different comparisons
cnt1 <- cnt1[c("X3.1.10.R1", "X3.1.10.R2", "X3.1.10.R3", "X3.3.1.R1.R3", "X3.3.1.R2", "X3.3.1.R3",
               "CR001.WT.R1", "CR001.WT.R2", "CR001.WT.R3", "CR007.R1", "CR007.R2", "CR007.R3")]

samp1 <- data.frame(sample = colnames(cnt1),
                          group = c(rep("3110_331", 6), rep("HD", 6)),
                          gender= c(rep("M",3), rep("F", 3), rep("M", 3), rep("F", 3)))

rownames(samp1) <- colnames(cnt1)


samp1$group <- as.factor(samp1$group)


dds <- DESeqDataSetFromMatrix(countData = cnt1,
                              colData = samp1,
                              design = ~ gender + group)

summary(cnt1)
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]
dds


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup = c( "group"), returnData = FALSE) #quality check


myalpha=.01

dds1 <- DESeq(dds)
res <- results(dds1, contrast = c('group', '3110_331', 'HD'))
summary(res)


sigCnt = sum(res$padj < myalpha, na.rm=TRUE)
cat("Count of differential atac peaks  ")
sigCnt
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < myalpha)

write.csv(as.data.frame(resSig), file="ipsc_111_131_versus_123_156_174_185_dp_results_2.csv")



pdf("rejections.pdf")
plot(metadata(res)$filterNumRej, type="b", ylab="Number of rejections", xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()

#MA plot
pdf("MAplot.pdf")
plotMA(res, ylim=c(-4,4))
abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
dev.off()

#Volcano plot?
pdf("Volcano.pdf")
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.01)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

#input dp_results
dp_results <- as.data.frame(resSig)

#switch to format input to Galaxy
library(tidyverse)
dp_results_new <- dp_results %>% 
  rownames_to_column(var = "row_name") %>% 
  separate(row_name,sep = ":",into = c("chrom", "Start", "End")) 

write.csv(dp_results_new, file="ipsc_111_131_versus_123_156_174_185_dp_results_new_2.csv", row.names=FALSE)
