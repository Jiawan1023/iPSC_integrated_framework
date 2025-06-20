library("tximport")
library("AnnotationHub")
library("ensembldb")
library( "biomaRt" )
library("dplyr")

dir <- "/input_path"
outdir <- "/output_path"
#example change sample name as needed 
fam1_111 <- list.files(path = dir, 
                   pattern= "fam1iPSC_111", 
                   full.names = FALSE)

fam1_156 <- list.files(path = dir, 
                       pattern= "fam1iPSC_156", 
                       full.names = FALSE)
#read in tx2gene R object previously saved
#tx2gene <- readRDS("/Users/jks6575/Dropbox/iPSC data summary/tx2gene.rds")


#if don't have tx2gene, need to generate.
tx2gene from AnnotationHub

hub = AnnotationHub()
edb <- query(hub, c("EnsDb", "sapiens", "109"))

hsapiens <- query(hub, c("EnsDb", "sapiens", "109"))[[1]]
keytypes(hsapiens)

txids <- keys(hsapiens, keytype = "TXID")
length(txids)

tx2gene <- select(hsapiens, txids, c("GENEID"), "TXID")

samples_all <- c(fam1_131, fam1_185) 
files <- file.path(dir, samples_all, "abundance.tsv") 
names(files) <- samples_all

#remove genes from  no chromosome from abundance files:


abundance.data <- lapply(files, read.table, header = TRUE)

mart <- useEnsembl(biomart = 'ensembl', 
                   dataset = 'hsapiens_gene_ensembl')

chrmap <- getBM( attributes = c("ensembl_transcript_id", "chromosome_name"),
                  filters = c("ensembl_transcript_id"),
                  values = gsub("\\..*", "", abundance.data[[1]]$target_id),
                  mart = mart )

values <- gsub("\\..*", "", abundance.data[[1]]$target_id)
nochrtxs <- values[!values %in% chrmap$ensembl_transcript_id]

#remove these from abundance table
abundance.data.nonochr <- abundance.data
for (i in 1:length(abundance.data.nonochr)) {
  abundance.data.nonochr[[i]] <- abundance.data.nonochr[[i]] %>%
    mutate(ensembl_transcript_id = gsub("\\..*", "", abundance.data.nonochr[[i]]$target_id))
  abundance.data.nonochr[[i]] <- subset(abundance.data.nonochr[[i]], !(abundance.data.nonochr[[i]]$ensembl_transcript_id) %in% nochrtxs)
}

#rerun getBM with corrected abundance table
chrmap <- getBM( attributes = c("ensembl_transcript_id", "chromosome_name"),
                 filters = c("ensembl_transcript_id"),
                 values = gsub("\\..*", "", abundance.data.nonochr[[1]]$target_id),
                 mart = mart )
dim(chrmap)[1] == dim(abundance.data.nonochr[[1]])[1]
#TRUE

#only keep normal chromosomes  and alterante contigs:

removethese <- unique(c(grep("H", chrmap$chromosome_name, value = TRUE), 
                        grep("K", chrmap$chromosome_name, value = TRUE), 
                        grep("G", chrmap$chromosome_name, value = TRUE), 
                        grep("MT", chrmap$chromosome_name, value = TRUE))) 
                    
#add a column with chromosome to abundance table for downstream filtering

abundance.data.withchr <- abundance.data.nonochr
for (i in 1:length(abundance.data.withchr)) {
  abundance.data.withchr[[i]] <- merge(abundance.data.withchr[[i]],
                                       chrmap,
                                       by = "ensembl_transcript_id",
                                       sort = FALSE)
}

abundance.data.nowchrom <- abundance.data.withchr
for (i in 1:length(abundance.data.nowchrom)) {
  abundance.data.nowchrom[[i]] <- abundance.data.nowchrom[[i]] %>% 
    dplyr::filter(!chromosome_name %in% removethese)
  write.table(abundance.data.nowchrom[[i]][,2:6],
              file.path(outdir, paste(names(files)[i], "_abundance.tsv", sep = "")),
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")
}

head(tx2gene)

tx2gene.nowchrom <- subset(tx2gene, tx2gene$TXID %in% abundance.data.nowchrom[[1]]$ensembl_transcript_id)

files.nowchrom <- list.files(path = outdir, 
                         pattern= "abundance.tsv", 
                         full.names = TRUE)

names4files <- list.files(path = outdir, 
                          pattern= "abundance.tsv", 
                          full.names = FALSE)
names4files <- gsub("_abundance.tsv", "", names4files)
names(files.nowchrom) <- names4files


txi <- tximport(files.nowchrom, type = "kallisto", tx2gene = tx2gene.nowchrom,
                ignoreTxVersion = TRUE)

#raw counts can be saved as dataframe using txi$counts
