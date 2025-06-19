#111 is p1C_007, 3110 is P2C_079 and 414 is P1C_077
#exmaple code

library(dplyr)

df1 <- fam1_111_polr3e_variants_fc0.5_density_input%>% mutate(Source = "111")
df2 <- fam3_3110_polr3e_variants_fc0.5_density_input%>% mutate(Source = "3110")
df3 <- fam4_414_polr3e_variants_fc0.5_density_input%>% mutate(Source = "414")

combined_df <- bind_rows(df1, df2) #change this for different comparisons 

library(ggplot2)

library(kSamples)

# Define custom colors for each source (new colors)
custom_colors <- c("111" = "#FF5733",   # Vibrant Red-Orange
                   "3110" = "#33FF57",  # Bright Green
                   "414" = "#3357FF")  # Deep Blue


# Define PDF file output with 3x3 inch dimensions
pdf("npc_111_3110_polr3e_density_plot.pdf", width = 3, height = 3)

# Create and save the density plot
ggplot(combined_df, aes(x = V3, fill = Source, color = Source)) +
  geom_density(alpha = 0.5) +  # Overlay density plots with transparency
  scale_fill_manual(values = custom_colors) +  # Apply custom fill colors
  scale_color_manual(values = custom_colors) +  # Apply custom line colors
  labs(title = "Density Plot of V3 for Different Sources", 
       x = "V3 Values", 
       y = "Density") +
  theme_minimal() +
  ylim(0, 0.125) + 
  xlim(0, 25)

# Close the PDF file
dev.off()

  # Anderson-Darling k-Sample Test
ad_test_result <- ad.test(list(df1$V3, df2$V3))
print(ad_test_result)
