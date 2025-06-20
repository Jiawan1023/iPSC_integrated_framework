library(IsoformSwitchAnalyzeR)

dir <- "C:/Users/Jiawan/Dropbox/Jiawan/Analysis/isoform/abundance_mn"

patterns <- c("111", "131", "185", "CR001WT", "CR007WT")

# List and concatenate files matching any pattern
samples_all <- unlist(lapply(patterns, function(p) list.files(path = dir, pattern = p, full.names = FALSE)))

# Generate full file paths and assign names
files <- setNames(file.path(dir, samples_all), samples_all)


kallistoQuant <- importIsoformExpression(sampleVector=files)

#condtion matrix
myDesign <- data.frame(
    sampleID = colnames(kallistoQuant$abundance)[-1],
    condition = c(rep("111_131_185", 8), rep("hd_control", 5))
)
myDesign

# Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix   = kallistoQuant$counts,
    isoformRepExpression = kallistoQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "C:/Users/Jiawan/Dropbox/Jiawan/Analysis/isoform/Homo_sapiens.GRCh38.112.chr.gtf",
    isoformNtFasta       = "C:/Users/Jiawan/Dropbox/Jiawan/Analysis/isoform/Homo_sapiens.GRCh38.cds.all.fa",
    ignoreAfterPeriod= TRUE,
    fixStringTieAnnotationProblem = TRUE,
    showProgress = FALSE
)

aSwitchList$isoformFeatures %>% head(2)


aSwitchList$orfAnalysis %>% head(2)

#prefilter
# geneExpressionCutoff = 1 
# isoformExpressionCutoff = 0 
# IFcutoff=0.01
# dIFcutoff = 0.1
# removeSingleIsoformGenes = TRUE
sar2 <- preFilter(aSwitchList)

# Perform test
sar3 <- isoformSwitchTestDEXSeq(sar2)

sar3$isoformFeatures %>% head(2)

sar3$isoformSwitchAnalysis %>% head()

setwd('C:/Users/Jiawan/Dropbox/Jiawan/Analysis/isoform')


output_prefix <- "mn_111_131_185_hd_control"
write.csv(sar3$isoformFeatures, file = paste0(output_prefix, "_isoformFeatures.csv"), row.names = FALSE)
write.csv(sar3$isoformSwitchAnalysis, file = paste0(output_prefix, "_isoformSwitchAnalysis.csv"), row.names = FALSE)
  
# Create the plot
p <- ggplot(data = sar3$isoformFeatures, aes(x = dIF, y = -log10(isoform_switch_q_value))) +
    geom_point(
    aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05), # default cutoff
    size = 1
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') + # default cutoff
    geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dashed') + # default cutoff
    scale_color_manual('Signficant\nIsoform Switch', values = c('black', 'red')) +
    labs(x = 'dIF', y = '-Log10 ( Isoform Switch Q Value )') +
    theme_bw()
  
# Save the plot as a PDF
pdf(paste0(output_prefix, "_volcano_plot.pdf"), height = 4, width = 5)
print(p) # Ensure the plot is printed to the PDF
dev.off()

#Predicting Alternative Splicing

sar4 <- analyzeAlternativeSplicing(sar3)
sar4$AlternativeSplicingAnalysis %>% head(2)


extractSplicingSummary(sar4)


q1 <- extractSplicingSummary(sar4)

q2 <- extractSplicingEnrichment(
    sar4,
    splicingToAnalyze='all',
    returnResult=FALSE,
    returnSummary=FALSE
)

q3 <- extractSplicingGenomeWide(
    sar4,
    splicingToAnalyze='all',
    returnResult=FALSE,
)


pdf(paste0(output_prefix, "_extractsplicingsummary_plot.pdf"), height = 4, width = 5)
print(q1) # Ensure the plot is printed to the PDF
dev.off()


pdf(paste0(output_prefix, "_extractsplicingenrichment_plot.pdf"), height = 4, width = 6)
    q2 # Ensure the plot is printed to the PDF
dev.off()

pdf(paste0(output_prefix, "_extractsplicinggenomwide_plot.pdf"), height = 4, width = 7)
    q3 # Ensure the plot is printed to the PDF
dev.off()
