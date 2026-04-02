
#Figure 1c
#raw data from Supplementary Table S3G and S3H
library(ComplexHeatmap)
library(ggplot2)

#input raw data as df

mat <- matrix(-log10(df$P.value), ncol = 1)
colnames(mat) <- "Motif Name"



ht <- Heatmap(
  mat,
  name = "-log10(p)",
  col = colorRamp2(c(2,5,10), c("white","pink","red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  width = unit(3, "cm")
)

#--------------------------------------------------------------------------------------
#Figure 4f
# Load the data
data <- read.csv("path")


data$line <- factor(data$line, levels = c("Order_of_lines")) # Define custom order of pathway descriptions


description_order <- c(
  "WP_IL18_SIGNALING_PATHWAY",
  "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING",
  "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
  "WP_INFLAMMATORY_RESPONSE_PATHWAY",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "KEGG_TGF_BETA_SIGNALING_PATHWAY",
  "REACTOME_SIGNALING_BY_NTRKS",
  "WP_PI3KAKT_SIGNALING_PATHWAY",
  "WP_FOCAL_ADHESION_PI3KAKTMTORSIGNALING_PATHWAY",
  "WP_WNT_SIGNALING",
  "REACTOME_SIGNALING_BY_WNT",
  "REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT",
  "WP_WNT_SIGNALING_WP428",
  "REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING",
  "REACTOME_FORMATION_OF_THE_BETA_CATENIN_TCF_TRANSACTIVATING_COMPLEX",
  "REACTOME_SIGNALING_BY_WNT_IN_CANCER",
  "WP_HIPPO_SIGNALING_REGULATION_PATHWAYS",
  "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
  "WP_GABA_RECEPTOR_SIGNALING",
  "WP_DNA_REPLICATION",
  "REACTOME_DNA_REPLICATION",
  "REACTOME_SYNTHESIS_OF_DNA",
  "NA"
)

# Set Description as a factor with the specified order
data$Description <- factor(data$Description, levels = description_order)


library(ggplot2)
library(RColorBrewer)
library(scales)

library(tidyr)
library(dplyr)

# Expand grid to ensure all tile positions exist
data_full <- expand_grid(
  line = unique(data$line),
  Description = unique(data$Description)
) %>%
  left_join(data, by = c("line", "Description"))  # join original NES values

# Create the palette function
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space = "Lab")

# Generate a continuous color scale from -3 to 3
color_scale <- myPalette(100)  # You can increase this for smoother gradients

p <- ggplot(data_full, aes(x = line, y = Description, fill = NES)) +
  geom_tile(color = "darkgrey", size = 0.5) +  # grey border for all tiles
  geom_text(aes(label = round(NES, 2)), color = "black", na.rm = TRUE) +
  scale_fill_gradientn(
    colors = myPalette(100),
    limits = c(-3, 3),
    na.value = "white"
  ) +
  labs(
    title = "Heatmap of NES by Cell Type",
    x = "Cell Type",
    y = "Pathway Description",
    fill = "NES"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave("heatmap_NES_new.pdf", plot = p, width = 12, height = 7)

#--------------------------------------------------------------------------------------
#Figure 4g


library(pheatmap)
pathway_data <- pathway_genes_heatmap_input #supplementary table S7F

tpm_data <- CRISPRa_all_counts_tpm


ensembl_ids <- pathway_data$ensembl_gene_id

filtered_tpm <- tpm_data[match(ensembl_ids, rownames(tpm_data)), , drop = FALSE]

# Convert to matrix
filtered_tpm_matrix <- as.matrix(filtered_tpm)

# Define column groups
column_groups <- factor(rep(c("npc_P1C_077_CRISPRa", "npc_P1C_007_CRISPRa", "npc_P2C_079_CRISPRa"), times = c(12, 12, 9)))

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

# Define the column groups 
column_groups <- factor(rep(c("npc_P1C_077_CRISPRa", "npc_P1C_007_CRISPRa", "npc_P2C_079_CRISPRa"), times = c(12, 12, 9)))


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
    Group = c("npc_P1C_077_CRISPRa" = "#E69F00", "npc_P1C_007_CRISPRa" = "#56B4E9", "npc_P2C_079_CRISPRa" = "#009E73"),
    Treatment = c("EV" = "#F8766D", "mosmo" = "#7CAE00", "uqcrc2" = "#00BFC4", "polr3e" = "#C77CFF")
)



color_palette <- colorRampPalette(c("navy", "white", "darkred"))(100)

gaps_col_positions <- c(12, 24)  # After column 12 and 24

# Plot heatmap
pheatmap(filtered_tpm_z_matrix, annotation_col = annotation_col, 
         annotation_colors = annotation_colors, scale = "none", cluster_rows = TRUE, 
         cluster_cols = FALSE, border_color = "black", gaps_col = gaps_col_positions, 
         color = color_palette,  clustering_method = "ward.D2")




#similar code was also used for Figure 6e, change the input



#--------------------------------------------------------------------------------------
#Supplementary Figure 3c

library(ggplot2)
df <- #input fisher's exact results" #Supplementary table S2E and S2F and S3C


df$Odds_ratio[is.infinite(df$Odds_ratio)] <- NA


df$Cell_type <- factor(df$Cell_type, levels = c("ipsc", "npc", "imn", "mn"))



# Plot with dark purple for NA (Inf) values
p <- ggplot(df, aes(x = Direction, y = Cell_type, fill = Odds_ratio)) +
    geom_tile(color = "black", size = 0.5) +  # Adds thicker black borders to the tiles
    scale_fill_gradient2(low = "navy", mid = "white", high = "darkred", midpoint = 1, limits = c(0, 6), na.value = "red4") +
    labs(x = "Cell Type", y = "Variant", fill = "Odds Ratio") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("CRISPR_deg_fisher_exact_heatmap.pdf", plot = p, width = 5, height = 5, units = "in")



#------------------------------------------------
#supplementary figure 5d
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







