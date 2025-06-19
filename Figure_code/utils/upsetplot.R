library(dplyr)
library(readr)
library(UpSetR)
library(GSEAB)

# Define the folder path where _down_to_up.csv files are stored
#input  folder_path

# List all _down_to_up.csv files in the directory
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty list to store ensembl_gene_id from each file
gene_lists <- list()

# Loop through each file to extract ensembl_gene_id
for (file in csv_files) {
  
  # Read the CSV file
  df <- read_csv(file)
  
  # Ensure 'ensembl_gene_id' column exists before processing
  if ("ensembl_gene_id" %in% colnames(df)) {
    
    # Extract unique ensembl_gene_id and store in list
    file_name <- tools::file_path_sans_ext(basename(file))  # Extract filename without extension
    gene_lists[[file_name]] <- unique(df$ensembl_gene_id)
    
  } else {
    message("Skipping file (No 'ensembl_gene_id' column found): ", file)
  }
}

# Convert list into a binary matrix for UpSet plot
gene_sets <- fromList(gene_lists)

# Generate the UpSet plot
upset(gene_sets, 
      nsets = length(gene_lists), 
      order.by = "freq", 
      mainbar.y.label = "Intersection Size", 
      sets.x.label = "Dataset Size")
