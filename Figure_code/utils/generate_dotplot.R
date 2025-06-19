#Figure 1e
# Load necessary libraries
library(tidyverse)
library(ggplot2)
expr_data <- npc_counts_tpm #input tpm

gene_row <- expr_data %>% filter(X == "ENSG00000204175") #input gene of interest


# Create long format data for plotting
plot_data <- data.frame(
  Expression = as.numeric(gene_row[1, 2:7]),  # Select expression values
  Group = factor(rep(c("Deletion", "Control"), each = 3), levels = c("Deletion", "Control"))
)

# Define custom colors to match your reference image
dot_colors <- c("Deletion" = "#BFCBE2", "Control" = "#2E297A")  # Light blue and dark blue

# Plot with transparent background, black axes, and custom dot colors
g <- ggplot(plot_data, aes(x = Group, y = Expression, color = Group)) +
  geom_jitter(width = 0.1, size = 3) +
  scale_color_manual(values = dot_colors) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 1.5, color = "black") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

# Display the plot
print(g)

#-------------------------------------------------------------
#fig 2b.c and suppfig 3c
library(ggplot2)
library(dplyr)

#input dataframe with tpm, condition_1 (del or nondel) and condition_2 (hit or nonhit).
#detailed information is provided in supplementary table 4
summary_data <- data %>%
  group_by(condition_1, condition_2) %>%
  summarise(
    mean_tpm = mean(tpm),
    se_tpm = sd(tpm) / sqrt(n())
  )

# Reorder condition_2 levels to display 'nonhit' first
summary_data$condition_2 <- factor(summary_data$condition_2, levels = c("nonhit", "hit"))

pdf(file = "NEK9_ipsc_plot.pdf", width = 3, height = 2)  # Save plot as PDF in landscape mode
ggplot(data, aes(x = factor(condition_2, levels = c("nonhit", "hit")), y = tpm, group = condition_1, color = condition_1, fill = condition_1)) +
  stat_summary(fun = median, geom = "line", aes(group = condition_1), size = 1) +  # Lines connect median values
  geom_jitter(shape = 21, size = 3, color = "black", width = 0.05) +  # Jittered data points to avoid overlap
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black", fatten = 3) +  # Median bar
  labs(
    title = "Interaction Plot: TPM by Condition with Median Bars",
    x = "Condition 2 (Nonhit / Hit)",
    y = "TPM",
    color = "Condition 1",
    fill = "Condition 1"
  ) +
  theme(
    plot.background = element_blank(),            # Remove outer background
    panel.background = element_blank(),           # Remove panel background
    panel.grid = element_blank(),                 # Remove grid lines
    axis.line = element_line(color = "black"),    # Add axis lines
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.key = element_blank()                  # Remove legend key background
  )

  dev.off()  # Close the PDF device


#------------------------------------------------------
#suppfig 1e, supplig 3b

# input compariosn name with desired order as ordered_comparison
#input publish dataset names as ordered_study_names

#input summary in Supplementary table S9F and Supplementary table S9G as df

df$test_file <- factor(df$test_file, levels = ordered_comparison)
df$reference_file   <- factor(df$reference_file, levels = ordered_study_names)
library(ggplot2)

bubble_plot <- ggplot(df, aes(x = reference_file, y = test_file)) +
    # Main bubbles
    geom_point(aes(size = overlap_size, fill = odds_ratio), shape = 21, color = "black") +
    
    # Bold red outline for Adjusted.P.value < 0.05
    geom_point(
        data = df[df$Adjusted.P.value < 0.05, ],
        aes(size = overlap_size),
        shape = 21,
        color = "darkred",
        stroke = 1.2
    ) +
    
    # Bold black outline for Adjusted.P.value > 0.05 but P.value< 0.05
    geom_point(
        data = df[df$Adjusted.P.value > 0.05 & df$P.value< 0.05, ],
        aes(size = overlap_size),
        shape = 21,
        color = "black",
        stroke = 1.2
    ) +
    
    # Color and size scales
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 1, name = "odds_ratio") +
    scale_size_continuous(range = c(2, 10), name = "Gene Number") +
    
    # Theme and axis labels
    theme(panel.grid.major = element_line(linetype = 1.0, color = "lightgrey"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank(),
          legend.position = "right",
          legend.box = "vertical") +
    ylab("Comparison") +
    xlab("Study Name")

bubble_plot
output_file <- "CRISPR_published_data_bubble_plot.pdf"  # Define the output file name
ggsave(output_file, plot = bubble_plot, width = 6, height = 5, units = "in")

#for disgenet
df$X <- factor(df$X, levels = ordered_comparison)

library(ggplot2)

bubble_plot <- ggplot(df, aes(x = Term, y = X)) +
    # Main bubbles
    geom_point(aes(size = Combined.Score, fill = Odds.Ratio), shape = 21, color = "black") +
    
    # Bold red outline for Adjusted.P.value < 0.05
    geom_point(
        data = df[df$Adjusted.P.value < 0.05, ],
        aes(size = Combined.Score),
        shape = 21,
        color = "darkred",
        stroke = 1.2
    ) +
    
    # Bold black outline for Adjusted.P.value > 0.05 but P.value< 0.05
    geom_point(
        data = df[df$Adjusted.P.value > 0.05 & df$P.value< 0.05, ],
        aes(size = Combined.Score),
        shape = 21,
        color = "black",
        stroke = 1.2
    ) +
    
    # Color and size scales
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 1, name = "odds_ratio") +
    scale_size_continuous(range = c(2, 10), name = "Combined.Score") +
    
    # Theme and axis labels
    theme(panel.grid.major = element_line(linetype = 1.0, color = "lightgrey"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank(),
          legend.position = "right",
          legend.box = "vertical") +
    ylab("Comparison") +
    xlab("Study Name")

bubble_plot
output_file <- "CRISPR_disgenet_bubble_plot.pdf"  # Define the output file name
ggsave(output_file, plot = bubble_plot, width = 6, height = 6, units = "in")

#---------------------------------------------------------------
#suppfig 8d and suppfig 11a
#input go data supplementary table 5 S5E and supplemntary table 8 S8F Select the representatoive ones
go_data$qvalue <- as.numeric(go_data$qvalue)
go_data$neg_log10_qvalue <- -log10(go_data$qvalue)


go_data$line <- factor(go_data$line, levels = c("111_mosmo", "111_polr3e", "111_uqcrc2", "3110_polr3e", "3110_uqcrc2", "414_mosmo", "414_polr3e", "414_uqcrc2"))

go_data$database_Description <- paste(go_data$database, go_data$Description, sep = "_")


# Create the bubble plot
ggplot(go_data, aes(x = line, y = database_Description, size = Count, fill = neg_log10_qvalue)) +
    geom_point(shape = 21, color = "black") +  # Shape 21 allows for a fill color
    scale_fill_gradient(low = "lightpink", high = "brown") +  # Customize color gradient for qvalue
    theme_minimal() +
    labs(
        x = "line",
        y = "Description",
        size = "Count",
        fill = "-log10(qvalue)",
        title = "GO Term Bubble Plot"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )














