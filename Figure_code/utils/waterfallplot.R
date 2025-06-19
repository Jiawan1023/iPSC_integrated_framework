#example waterfall plot code
#core code is same for every waterfall plots in the paper. 

library(ggplot2)
library(dplyr)
library(ggrepel)

# Load data
gsea_data <- CRISPR_GSEA_summary

# Order data by Count or NES
gsea_data <- gsea_data %>% arrange(desc(NES))
gsea_data$Rank <- 1:nrow(gsea_data)  # Create a rank for x-axis

# Define shape mapping for each cell type
shape_map <- c(
  "ipsc" = 16,  # dot
  "npc"  = 17,  # triangle
  "imn"  = 15,  # square
  "mn"   = 18   # diamond
)

gsea_data <- gsea_data %>%
  mutate(
    color_custom = case_when(
      is.na(feature) ~ "grey70",
      feature == "signaling_pathway" & cell_type == "ipsc" ~ "#fc8d62",
      feature == "signaling_pathway" & cell_type == "npc"  ~"#8da0cb",
      feature == "signaling_pathway" & cell_type == "imn"  ~ "#e78ac3",
      feature == "signaling_pathway" & cell_type == "mn"   ~ "#66c2a5",
      TRUE ~ "black"  # fallback color
    )
  )

# Generate plot
gsea_plot <- ggplot(gsea_data, aes(x = Rank, y = NES)) +
  geom_point(aes(color = color_custom, shape = cell_type), size = 1.5, alpha = 1) +
  scale_color_identity() +
  scale_shape_manual(values = shape_map) +
  geom_text_repel(aes(label = ifelse(
    Description %in% c(
    "WP_GABA_RECEPTOR_SIGNALING", "REACTOME_GABA_SYNTHESIS_RELEASE_REUPTAKE_AND_DEGRADATION", "WP_PI3KAKT_SIGNALING_PATHWAY",
      "KEGG_TGF_BETA_SIGNALING_PATHWAY"
    ),
    Description, ""
  )), size = 3.5, max.overlaps = 60) +
  theme_minimal() +
  labs(x = "Rank", y = "NES", title = "GO Enrichment") +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

gsea_plot
