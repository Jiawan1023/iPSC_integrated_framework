#Figure 1c
# input supplementary table S2J as df

library(ggplot2)
df$Odds_ratio[is.infinite(df$Odds_ratio)] <- NA


df$Cell_type <- factor(df$Cell_type, levels = c("ipsc", "npc", "imn", "mn"))
df$Variant <- factor (df$Variant, levels = rev(c("snv_coding", "snv_noncoding", "STR_coding", "STR_noncoding", "deletion_coding", "deletion_noncoding", "duplication_coding", "duplication_noncoding")))


# Plot 
p <- ggplot(df, aes(x = Cell_type, y = Variant, fill = Odds_ratio)) +
    geom_tile(color = "black", size = 0.5) +  # Adds thicker black borders to the tiles
    scale_fill_gradient2(low = "navy", mid = "white", high = "darkred", midpoint = 1, limits = c(0, 6), na.value = "red4") +
    labs(x = "Cell Type", y = "Variant", fill = "Odds Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#input supplementary table S2K as df
df$cell_type <- factor(df$cell_type, levels = c("ipsc", "npc", "imn", "mn"))

df$Variant <- factor(
  df$Variant,
  levels = rev(c("splice", "intron", "STR_intronic", "deletion_intron", "duplication_intron"))
)

# Plot with reordered y-axis
p <- ggplot(df, aes(x = cell_type, y = Variant, fill = Odds_ratio)) + 
  geom_tile(color = "black", size = 0.5) +  # Adds thicker black borders to the tiles
  scale_fill_gradient2(
    low = "navy", mid = "white", high = "darkred", 
    midpoint = 1, limits = c(0, 8), na.value = "red4"
  ) + 
  labs(x = "Cell Type", y = "Variant", fill = "Odds Ratio") + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#------------------------------------------------------

#Figure 1d
#input supplementary table S3I as df
df$Odds_ratio[is.infinite(df$Odds_ratio)] <- NA


df$cell_type <- factor(df$cell_type, levels = c("ipsc", "npc"))
df$Variant <- factor (df$Variant, levels = rev(c("utr5", "utr3", "intron", "upstream", "downstream",
                                                 "STR_UTR5", "STR_UTR3", "STR_intronic", "STR_upstream", "STR_downstream",
                                                 "deletion_utr5", "deletion_utr3", "deletion_intron", "deletion_upstream", "deletion_downstream",
                                                 "duplication_utr5", "duplication_utr3", "duplication_intron","duplication_upstream", "duplication_downstream")))


# Plot with dark purple for NA (Inf) values
p <- ggplot(df, aes(x = cell_type, y = Variant, fill = Odds_ratio)) +
    geom_tile(color = "black", size = 0.5) +  # Adds thicker black borders to the tiles
    scale_fill_gradient2(low = "navy", mid = "white", high = "darkred", midpoint = 1, limits = c(0, 12), na.value = "red4") +
    labs(x = "Cell Type", y = "Variant", fill = "Odds Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Suppfig 1g ans Suppfig 11a
#similar code, source data is Supplementary Table 3 S3D for suppfig 1g
#input Supplementary Table 10 S8E as df for suppfig 11a


#--------------------------------------------------------------------------------------
#Figure 4g and figure 6e
library(pheatmap)
pathway_data <- pathway_genes_heatmap_input #supplementary table S5G
#change the pathway data to TF candidate genes for figure 6e.
tpm_data <- CRISPRa_all_counts_tpm


ensembl_ids <- pathway_data$ensembl_gene_id

filtered_tpm <- tpm_data[match(ensembl_ids, rownames(tpm_data)), , drop = FALSE]

# Convert to matrix
filtered_tpm_matrix <- as.matrix(filtered_tpm)

# Convert to matrix
filtered_tpm_matrix <- as.matrix(filtered_tpm)

# Define column groups
column_groups <- factor(rep(c("npc_111_CRISPRa", "npc_414_CRISPRa", "npc_3110_CRISPRa"), times = c(12, 12, 9)))

# Function to normalize by group mean and then scale within each group
normalize_and_scale <- function(mat, groups) {
  unique_groups <- unique(groups)
  normalized_mat <- mat  # Placeholder for results
  
  for (group in unique_groups) {
    cols <- which(groups == group)  # Get column indices for the group
    group_data <- mat[, cols, drop = FALSE]  # Extract relevant columns
    
    # Normalize by group mean (divide each row by mean of group)
    group_means <- rowMeans(group_data, na.rm = TRUE)
    group_data <- sweep(group_data, 1, group_means, "/")  # Divide each row by its mean
    
    # Scale each row within the group
    group_data <- t(scale(t(group_data)))  # Z-score per row within group
    
    # Store back results
    normalized_mat[, cols] <- group_data
  }
  
  return(normalized_mat)
}

# Apply normalization and scaling
filtered_tpm_z_matrix <- normalize_and_scale(filtered_tpm_matrix, column_groups)

# Convert back to dataframe (optional)
filtered_tpm_z <- as.data.frame(filtered_tpm_z_matrix)

# Ensure only matching rows are used
filtered_tpm_z <- filtered_tpm_z[rownames(filtered_tpm_z) %in% pathway_data$ensembl_gene_id, ]

# Replace Ensembl IDs with HGNC symbols
rownames(filtered_tpm_z) <- pathway_data$hgnc_symbol[match(rownames(filtered_tpm_z), pathway_data$ensembl_gene_id)]

# Convert to matrix for heatmap plotting
filtered_tpm_z_matrix <- as.matrix(filtered_tpm_z)

# Create a column annotation for grouping
annotation_col <- data.frame(Group = column_groups)
rownames(annotation_col) <- colnames(filtered_tpm_z_matrix)

# Define the column groups (pre-existing)
column_groups <- factor(rep(c("npc_111_CRISPRa", "npc_414_CRISPRa", "npc_3110_CRISPRa"), times = c(12, 12, 9)))
#111 is P1C_077, 414 is P1C_007, 3110 is P2C_079.

# Define sample annotations
sample_annotation <- c(
    rep("EV", 3), rep("mosmo", 3), rep("uqcrc2", 3), rep("polr3e", 3),  # npc_111_CRISPRa (12)
    rep("EV", 3), rep("mosmo", 3), rep("uqcrc2", 3), rep("polr3e", 3),  # npc_414_CRISPRa (12)
    rep("EV", 3), rep("uqcrc2", 3), rep("polr3e", 3)  # npc_3110_CRISPRa (9)
)

# Create a dataframe for annotations
annotation_col <- data.frame(
    Group = column_groups,
    Treatment = sample_annotation
)
rownames(annotation_col) <- colnames(filtered_tpm_z_matrix)

# Define colors for annotations
annotation_colors <- list(
    Group = c("npc_111_CRISPRa" = "#E69F00", "npc_414_CRISPRa" = "#56B4E9", "npc_3110_CRISPRa" = "#009E73"),
    Treatment = c("EV" = "#F8766D", "mosmo" = "#7CAE00", "uqcrc2" = "#00BFC4", "polr3e" = "#C77CFF")
)



color_palette <- colorRampPalette(c("navy", "white", "darkred"))(100)

gaps_col_positions <- c(12, 24)  # After column 12 and 24

# Plot heatmap
pheatmap(filtered_tpm_z_matrix, annotation_col = annotation_col, 
         annotation_colors = annotation_colors, scale = "none", cluster_rows = TRUE, 
         cluster_cols = FALSE, border_color = "black", gaps_col = gaps_col_positions, 
         color = color_palette,  clustering_method = "ward.D2")


#------------------------------------------------
#supplementary figure 4c
#input tpm counts for heatmap

                #nkx2-1,          foxg1,              dlx1,             dlx2,                dlx5,            dlx6,             lhx6,               lhx8,              gad2,            sox6
gaba_gene <- c("ENSG00000136352", "ENSG00000176165", "ENSG00000144355" ,"ENSG00000115844", "ENSG00000105880","ENSG00000006377", "ENSG00000106852", "ENSG00000162624", "ENSG00000136750","ENSG00000110693")

df<- counts_for_heatmap

df <- tibble::rownames_to_column(df, var = "Gene")


#extract gaba_gene_row
library(dplyr)
df <- df %>%
  filter(Gene %in% gaba_gene)

#rowname counts

# Convert 'Gene' column to row names
row.names(df) <- df$Gene

# Remove the 'Gene' column from the data frame
df <- df[, -which(names(df) == "Gene")]


write.csv(df,
          file = file.path(outdir, "counts_for_heatmap_gaba_gene.csv"),
          quote = FALSE,
          row.names = TRUE,
          col.names = TRUE)

library(pheatmap)

annotation_col <- data.frame(Condition = c(rep("ipsc",18), rep("npc",18), rep("imn", 18), rep("mn",17)))

rownames(annotation_col) <- colnames(counts_for_heatmap_gaba_gene)

annotation_col$Family <- c(rep("CRISPR_deletion",3), rep("family4",6), rep("controls", 9), rep("CRISPR_deletion",3), rep("family4",6), rep("controls", 9),rep("CRISPR_deletion",3), rep("family4",6), rep("controls", 9),rep("CRISPR_deletion",3), rep("family4",6), rep("controls", 8))
p <- pheatmap(df, 
              cluster_rows = FALSE,                       
              cluster_cols = FALSE, 
              show_rownames = TRUE, 
              show_colnames = TRUE,
              annotation_col = annotation_col)

#family4 is GL_007.
#controls contains CRISPR control and HD_01 and HD_02








