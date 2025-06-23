library(enrichR)
setEnrichrSite("Enrichr")

# Set up the database
dbs <- c("DisGeNET")

# Specify the folder containing the files
input_folder <- "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_enrichment/DEGs_for_fisher/test_gene/fam3"  # Replace with your folder path
output_folder <- "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_enrichment/disgenet_enrichment/result"  # Replace with your output folder path

csv_files <- list.files(input_folder, full.names = TRUE, pattern = "//.csv$")  # Adjust pattern for file types
for (csv_file in csv_files) {
  # Read the input file
  test_gene <- read.csv(csv_file)

  # Add a significance column based on criteria
  test_gene[which(test_gene$log2FoldChange >= 1 & test_gene$padj < 0.05), 'sig'] <- 'high'
  test_gene[which(test_gene$log2FoldChange <= -1 & test_gene$padj < 0.05), 'sig'] <- 'high'
  test_gene[which(abs(test_gene$log2FoldChange) <= 1 | test_gene$padj >= 0.05), 'sig'] <- 'none'

  # Filter for significant genes
  test_gene_select <- subset(test_gene, sig %in% c('high'))

  # Run Enrichr analysis
  enriched <- enrichr(test_gene_select$external_gene_name, dbs)

  # Convert the results to a data frame
  df <- as.data.frame(enriched[["DisGeNET"]])

  # Generate output file name based on input file name
  input_file_name <- basename(csv_file)  # Extract input file name
  output_file_name <- sub("//.csv$", "_enrichment.csv", input_file_name)  # Append '_enrichment'
  output_file <- file.path(output_folder, output_file_name)  # Full output path

  # Write results to a CSV file
  write.csv(df, output_file, row.names = FALSE)
}


# Define the input folder path
input_folder <- "C:/Users/Jiawan/Dropbox/Jiawan/Analysis/Genetic_enrichment/disgenet_enrichment/result/archive"

# Specify the terms to filter
filter_terms <- c("Schizophrenia", "Bipolar Disorder", "Depressive disorder", 
                  "Epileptic Seizures", "Autism Spectrum Disorders", "Anxiety Disorders")

# Initialize an empty list to store filtered data
filtered_list <- list()

# List all CSV files in the folder
csv_files <- list.files(input_folder, full.names = TRUE, pattern = "\\.csv$")

# Loop through each CSV file
for (csv_file in csv_files) {
    # Extract the file name (without extension)
    file_name <- tools::file_path_sans_ext(basename(csv_file))
    
    # Load the data frame
    df <- read.csv(csv_file)
    
    # Check if the Term column exists
    if (!"Term" %in% colnames(df)) {
        cat("Skipping file (no Term column):", csv_file, "\n")
        next
    }
    
    # Filter rows where the Term column matches the filter terms
    filtered_df <- subset(df, Term %in% filter_terms)
    
    # Add the filtered data frame to the list if it has any rows
    if (nrow(filtered_df) > 0) {
        filtered_list[[file_name]] <- filtered_df
    }
}

# Combine all filtered data frames into one
combined_filtered_df <- do.call(rbind, filtered_list)

# View the resulting data frame
print(head(combined_filtered_df))

# Save the combined filtered data frame to a CSV file
output_file <- "combined_results_disgenet.csv"
write.csv(combined_filtered_df, output_file, row.names = TRUE)

cat("Filtered results saved to:", output_file, "\n")
