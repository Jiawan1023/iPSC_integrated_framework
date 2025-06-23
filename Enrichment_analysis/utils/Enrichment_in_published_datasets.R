input_folder <- "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_enrichment/DEGs_for_fisher/test_gene/CRISPR/"

reference_folder <- "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_enrichment/DEGs_for_fisher/reference_gene/"




# List all CSV files in the input folder
csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

# List all Excel files in the reference folder
xlsx_files <- list.files(reference_folder, pattern = "\\.xlsx$", full.names = TRUE)

# Initialize an empty dataframe to store results
results_df <- data.frame(
  test_file = character(),
  reference_file = character(),
  p_value = numeric(),
  odds_ratio = numeric(),
  overlap_size = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each Excel file for reference genes
for (xlsx_file in xlsx_files) {
  # Read the Excel file
  reference_gene <- read_excel(xlsx_file)
  
  # Check if the reference_gene file contains the expected column
  if (!"gene_symbol" %in% colnames(reference_gene)) {
    print(paste("Skipping file:", xlsx_file, "- Missing 'gene_symbol' column"))
    next
  }
  
  # Loop through each CSV file for test genes
  for (csv_file in csv_files) {
    # Read the CSV file
    test_gene <- read.csv(csv_file)
    
    # Process the data
    test_gene[which(test_gene$log2FoldChange >= 1 & test_gene$padj < 0.05), 'sig'] <- 'high'
    test_gene[which(test_gene$log2FoldChange <= -1 & test_gene$padj < 0.05), 'sig'] <- 'high'
    test_gene[which(abs(test_gene$log2FoldChange) <= 1 | test_gene$padj >= 0.05), 'sig'] <- 'none'
    
    test_gene_select <- subset(test_gene, sig %in% c('high'))
    
    # Set population gene
    custom_universe <- unique(test_gene$external_gene_name)
    
    # Perform Fisher's exact test
    genome_size <- length(custom_universe)
    go.obj <- newGeneOverlap(unique(test_gene_select$external_gene_name), 
                             reference_gene$gene_symbol,
                             genome.size = genome_size)
    go.obj <- testGeneOverlap(go.obj)
    
    # Extract results from go.obj
    p_value <- go.obj@pval
    odds_ratio <- go.obj@odds.ratio
    overlap_size <- length(go.obj@intersection)
    
    # Append results to the dataframe
    results_df <- rbind(results_df, data.frame(
      test_file = basename(csv_file),
      reference_file = basename(xlsx_file),
      p_value = p_value,
      odds_ratio = odds_ratio,
      overlap_size = overlap_size,
      stringsAsFactors = FALSE
    ))
  }
}

# Add FDR_BH column to the results
results_df$FDR_BH <- p.adjust(results_df$p_value, method = "BH")

# Print the results dataframe
print(results_df)

# Optionally save the results to a CSV file
write.csv(results_df, file = "CRISPR_fisher_summary.csv", row.names = FALSE)
