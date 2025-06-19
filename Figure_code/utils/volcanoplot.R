library(ggplot2)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)
library(patchwork)
library(ggrastr)
# Genes to label
genes_to_label <- c('MOSMO', 'POLR3E', 'CDR2', 'UQCRC2', 'PDZD9', 'EEF2K', 'VWA3A')
#input deseq2 files (all in supplementary tables) as df.
# Define thresholds
logfc_threshold <- 0.5
pval_threshold <- 0.05

# Add significance column
df$significance <- "Not Significant"
df$significance[df$padj < pval_threshold & df$log2FoldChange > logfc_threshold] <- "Upregulated"
df$significance[df$padj < pval_threshold & df$log2FoldChange < -logfc_threshold] <- "Downregulated"

# Subset for labeled genes
df_label <- subset(df, external_gene_name %in% genes_to_label)

# Volcano plot with clear labels
p<- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
 geom_point_rast(alpha = 0.6) +
  # Highlight labeled points with bigger dots and black border
  geom_point(data = df_label, aes(x = log2FoldChange, y = -log10(padj)),
             color = "black", fill = "black", shape = 21, size = 1, stroke = 1) +
  # Add text labels
  geom_text_repel(data = df_label, aes(label = external_gene_name),
                  size = 4, fontface = "bold", box.padding = 0.5, point.padding = 0.4,
                  max.overlaps = 100) +
  scale_color_manual(values = c("Upregulated" = "darkred", "Downregulated" = "navyblue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot with Highlighted Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed")
ggsave("mn_CRISPR_volcano_plot_labeled.pdf", plot = p, width = 5, height = 5)


