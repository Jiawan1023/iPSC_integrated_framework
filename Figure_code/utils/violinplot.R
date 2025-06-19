#example code. core codes for every violin plots in the paper is same.

df.summary <- df %>%
     group_by(condition) %>% 
     summarise(
         sd = sd(Sox2, na.rm = TRUE),
         Sox2 = mean(Sox2)
     )
 df.summary


library(ggplot2)
library(ggpubr)
library(RColorBrewer)

p <- ggplot(df, aes(condition, Sox2, fill = condition)) + 
    geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +  # Box plot inside violin
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) +  # Custom fill colors
    scale_color_manual(values = c("#00AFBB", "#E7B800")) + 
    theme(
        legend.position = "top",
        panel.background = element_rect(fill = "white"),  # White background for the panel
        plot.background = element_rect(fill = "white")  # White background for the plot
    ) +
    scale_fill_brewer(palette = "Set2") +
    stat_compare_means(method = "t.test", label = "p.format", vjust = -6)  # Add t-test comparison

ggsave("Sox2_violin_plot.svg", plot = p, width = 6, height = 5, dpi = 300, device = "svg")
