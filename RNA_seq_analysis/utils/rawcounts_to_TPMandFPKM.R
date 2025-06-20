#count gene exon length
#for isoform, use the right reference file.
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("gencode.v42.chr_patch_hapl_scaff.basic.annotation.gtf",format="gtf")
 exons_gene <- exonsBy(txdb, by ="gene")
 exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
 exons_gene_lens1 <- as.data.frame(exons_gene_lens)
 exons_gene_lens1 <- t(exons_gene_lens1)
write.csv(exon_gene_lens1,file = "exon_gene_lens.csv",row.names = T)

#import raw counts and exon gene lens read.csv("")
merged_dataset <- merge(rawcounts, exons_gene_lens)
write.csv(merged_dataset,file = "merge_rawcounts.csv",sep = "\t")
mycounts<-read.csv("merge_rawcounts.csv")
rownames(mycounts)<-mycounts[,2] #use ensemble ids as rownames
mycounts<-mycounts[,-1] #delete unnecessary rows
names(mycounts)[8] <- "Length" #name the gene lenghth column 
kb <- mycounts$Length / 1000
countdata <- mycounts[,2:7] #only extract raw counts data
rpk <- countdata / kb
tpm <- t(t(rpk)*1e6 / colSums(rpk))
write.csv(tpm,file = "rawcounts_tpm.csv",row.names = T)
fpkm <- t(t(rpk)/colSums(countdata)*1e6) 
write.csv(fpkm,file = "rawcounts_fpkm.csv",row.names = T)
