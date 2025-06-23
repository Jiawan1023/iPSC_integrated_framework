#for gene
import pandas as pd
import glob
import os
import re
#input atacannotation files
# Define directories
csv_path = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/result/ipsc/gene/atac"  # Change to your CSV folder
annotated_path = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/result/ipsc/gene/annotated"  # Change to your Annotated folder

# Get all CSV files
csv_files = glob.glob(os.path.join(csv_path, "*_annotation.csv"))
print(csv_files)

# Get all Annotated files
annotated_files = glob.glob(os.path.join(annotated_path, "annotated_*.csv"))


# Function to extract key from CSV filename
def extract_csv_key(filename):
    match = re.search(r"ipsc_(.*?)_versus", filename)
    return match.group(1) if match else None

# Function to extract key from Annotated filename
def extract_annotated_key(filename):
    match = re.search(r"annotated_significant_result_(.*?)_versus", filename)
    return match.group(1) if match else None

# Create dictionaries to store data
csv_data = {}
csv_filenames = {}
for file in csv_files:
    key = extract_csv_key(os.path.basename(file))
    if key:
        csv_data[key] = pd.read_csv(file)
        csv_filenames[key] = os.path.basename(file)  # Store filename

annotated_data = {}
annotated_filenames = {}
for file in annotated_files:
    key = extract_annotated_key(os.path.basename(file))
    if key:
        if key not in annotated_data:
            annotated_data[key] = []  # Allow multiple files per key
            annotated_filenames[key] = []  # Store multiple filenames
        annotated_data[key].append(pd.read_csv(file))
        annotated_filenames[key].append(os.path.basename(file))



# Pair CSV and Annotated files based on extracted keys
paired_files = {}
print("\n**Paired Files:**\n")
output_path = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/result/ipsc/gene/atac_annotated"
os.makedirs(output_path, exist_ok=True)  # Ensure the output directory exists
for key in csv_data.keys():
    if key in annotated_data:
        df_csv = csv_data[key]
        csv_file_name = csv_filenames[key]
        
        paired_files[key] = []  # Allow multiple annotated files per key

        # Iterate over all matching annotated files
        for idx, df_annotated in enumerate(annotated_data[key]):
            annotated_file_name = annotated_filenames[key][idx]

            # Ensure 'ensembl_gene_id' exists in both DataFrames
            if 'ensembl_gene_id' in df_csv.columns and 'ensembl_gene_id' in df_annotated.columns:
                # Create a set of Ensembl IDs from the CSV file
                ensembl_set = set(df_csv['ensembl_gene_id'])

                # Add atac_genename column (1 if match, 0 otherwise)
                df_annotated['atac_genename'] = df_annotated['ensembl_gene_id'].apply(lambda x: "1" if x in ensembl_set else "0")

                # Move atac_genename to the first column
                cols = ['atac_genename'] + [col for col in df_annotated.columns if col != 'atac_genename']
                df_annotated = df_annotated[cols]

                # Store the updated DataFrame
                paired_files[key].append(df_annotated)
                # Save the updated annotated file
                output_filename = f"atac_{annotated_file_name}"
                output_filepath = os.path.join(output_path, output_filename)
                df_annotated.to_csv(output_filepath, index=False)


                # Print the file pair
                print(f"CSV File: {csv_file_name}  <-->  Annotated File: {annotated_file_name}")

print("\nTotal unique CSV keys matched:", len(paired_files))



variant_path = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/result/ipsc/gene/variant"

# Get all Variant files starting with "variant_"
variant_files = glob.glob(os.path.join(variant_path, "filtered_variant_*.csv"))

# Function to process variant file
def process_variant_file(file_path):
    df = pd.read_csv(file_path)

    def extract_vid(value):
        """Extracts chromosome (first part before '_')"""
        return value.split('_')[0]

    def extract_start(value):
        """Extracts first number after '_'"""
        match = re.search(r'_(\d+)', value)
        return int(match.group(1)) if match else None

    def extract_end(value):
        """Extracts second number after second '_', otherwise uses start"""
        match = re.search(r'_(\d+)_(\d+)', value)
        if match:
            return int(match.group(2))
        start_match = re.search(r'_(\d+)', value)
        return int(start_match.group(1)) if start_match else None

    # Apply transformations
    df["chrom"] = df["vid"].apply(extract_vid)
    df["Start"] = df["vid"].apply(extract_start)
    df["End"] = df["vid"].apply(extract_end)

    # Save the modified file
    output_file = os.path.join(variant_path, "modified_" + os.path.basename(file_path))
    df.to_csv(output_file, index=False)
    print(f"Processed and saved: {output_file}")

# Process each variant file
for variant_file in variant_files:
    process_variant_file(variant_file)



csv_results_path = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/result/ipsc/gene/atac"
variant_path = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/result/ipsc/gene/variant"
atac_annotated_path = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/result/ipsc/gene/atac_annotated"

# Get all CSV files ending with "_results_new.csv"
csv_results_files = glob.glob(os.path.join(csv_results_path, "*_results_new.csv"))

# Get all Annotated files starting with "annotated_significant_result"
atac_annotated_files = glob.glob(os.path.join(atac_annotated_path, "atac_annotated_significant_result_*.csv"))

# Get all Variant files starting with "variant_"
variant_files = glob.glob(os.path.join(variant_path, "modified_filtered_variant_*.csv"))

# Function to extract key from results_new CSV filename
def extract_results_key(filename):
    match = re.search(r"ipsc_(.*?)_versus", filename)
    return match.group(1) if match else None

# Function to extract key from Annotated filename
def extract_annotated_key(filename):
    match = re.search(r"atac_annotated_significant_result_(.*?)_versus", filename)
    return match.group(1) if match else None

# Function to extract key from Variant filename
def extract_variant_key(filename):
    match = re.search(r"filtered_variant_(.*?)_versus", filename)
    return match.group(1) if match else None

# Create dictionaries to store filenames
csv_results_filenames = {}
for file in csv_results_files:
    key = extract_results_key(os.path.basename(file))
    if key:
        csv_results_filenames[key] = os.path.basename(file)

atac_annotated_filenames = {}
for file in atac_annotated_files:
    key = extract_annotated_key(os.path.basename(file))
    if key:
        if key not in atac_annotated_filenames:
            atac_annotated_filenames[key] = []
        atac_annotated_filenames[key].append(os.path.basename(file))

variant_filenames = {}
for file in variant_files:
    key = extract_variant_key(os.path.basename(file))
    if key:
        variant_filenames[key] = os.path.basename(file)

# Pair Variant files with Results and Annotated files
print("\n**Paired Files (variant, results_new.csv & Annotated files):**\n")

paired_files = {}
for key in csv_results_filenames.keys():
    if key in atac_annotated_filenames and key in variant_filenames:
        csv_results_file_name = csv_results_filenames[key]
        variant_file_name = variant_filenames[key]

        paired_files[key] = atac_annotated_filenames[key]  # Store multiple annotated files for the same key

        for atac_annotated_file_name in atac_annotated_filenames[key]:
            print(f"Variant File: {variant_file_name}  <-->  Results CSV File: {csv_results_file_name}  <-->  Annotated File: {atac_annotated_file_name}")

print("\nTotal unique CSV keys matched:", len(paired_files))


import os
import re
import glob
import pandas as pd
import pybedtools

# Get all CSV files
csv_results_files = glob.glob(os.path.join(csv_results_path, "*_results_new.csv"))
variant_files = glob.glob(os.path.join(variant_path, "modified_filtered_variant_*.csv"))
atac_annotated_files = glob.glob(os.path.join(atac_annotated_path, "atac_annotated_significant_result_*.csv"))

# Function to extract key from filenames
def extract_key(filename, pattern):
    match = re.search(pattern, filename)
    return match.group(1) if match else None

# Create dictionaries to store filenames grouped by key
csv_results_dict = {extract_key(f, r"ipsc_(.*?)_versus"): f for f in csv_results_files if extract_key(f, r"ipsc_(.*?)_versus")}
variant_dict = {}
atac_annotated_dict = {}

# Group variant_files by key
for f in variant_files:
    key = extract_key(f, r"modified_filtered_variant_(.*?)_versus")
    if key:
        if key not in variant_dict:
            variant_dict[key] = []
        variant_dict[key].append(f)

# Group atac_annotated_files by key
for f in atac_annotated_files:
    key = extract_key(f, r"atac_annotated_significant_result_(.*?)_versus")
    if key:
        if key not in atac_annotated_dict:
            atac_annotated_dict[key] = []
        atac_annotated_dict[key].append(f)

# Debug: Print detected files
print("\nDetected Files:")
print(f"Results: {csv_results_dict}")
print(f"Variants: {variant_dict}")
print(f"Annotated: {atac_annotated_dict}")

# Process each csv_results_file and pair with matching variant_files and atac_annotated_files
for results_key, results_file in csv_results_dict.items():
    if results_key not in variant_dict or results_key not in atac_annotated_dict:
        print(f"Skipping {results_key} - Missing variant or annotated files.")
        continue

    variant_files_for_key = variant_dict[results_key]
    atac_annotated_files_for_key = atac_annotated_dict[results_key]

    print(f"\nProcessing: Results ({results_file})")
    print(f"Paired Variants: {variant_files_for_key}")
    print(f"Paired Annotated: {atac_annotated_files_for_key}")

    # Load results file
    try:
        results_df = pd.read_csv(results_file)
    except Exception as e:
        print(f"Error loading results file {results_file}: {e}")
        continue

    # Ensure required columns exist in results file
    required_cols_results = {"chrom", "Start", "End"}
    if not required_cols_results.issubset(results_df.columns):
        print(f"Skipping {results_key} - Missing required columns in results file.")
        continue

    # Convert results to BED format for pybedtools
    results_bed = results_df[["chrom", "Start", "End"]].dropna()
    results_bed = results_bed.astype({"chrom": str, "Start": int, "End": int})
    #results_bed["Start"] = results_bed["Start"] - 500  # Extend start by 500 bp
    #results_bed["End"] = results_bed["End"] + 500  # Extend end by 500 bp
    results_bedtool = pybedtools.BedTool.from_dataframe(results_bed)

    # Process each paired variant file
    for variant_file in variant_files_for_key:
        try:
            variant_df = pd.read_csv(variant_file)
        except Exception as e:
            print(f"Error loading variant file {variant_file}: {e}")
            continue

        # Ensure required columns exist in variant file
        required_cols_variant = {"chrom", "Start", "End", "gene_id"}
        if not required_cols_variant.issubset(variant_df.columns):
            print(f"Skipping {variant_file} - Missing required columns in variant file.")
            continue

        # Convert variant to BED format for pybedtools
        variant_bed = variant_df[["chrom", "Start", "End"]].dropna()
        variant_bed = variant_bed.astype({"chrom": str, "Start": int, "End": int})
        variant_bedtool = pybedtools.BedTool.from_dataframe(variant_bed)

        # Perform intersection
        intersected = variant_bedtool.intersect(results_bedtool, wa=True)

        # Print intersected regions
        print(f"\nIntersected regions between {results_file} and {variant_file}:")
        for line in intersected:
            print(f"Chrom: {line.chrom}, Start: {line.start}, End: {line.end}")

        # Extract intersected regions and corresponding gene_id
        intersected_regions = {}
        for line in intersected:
            chrom, start, end = line.chrom, str(line.start), str(line.end)
            # Find the gene_id for the intersected region
            gene_id = variant_df[
                (variant_df["chrom"] == chrom) &
                (variant_df["Start"] == int(start)) &
                (variant_df["End"] == int(end))
            ]["gene_id"].values
            if len(gene_id) > 0:
                intersected_regions[(chrom, start, end)] = gene_id[0]

        # Process each paired annotated file
        for atac_annotated_file in atac_annotated_files_for_key:
            try:
                atac_annotated_df = pd.read_csv(atac_annotated_file)
            except Exception as e:
                print(f"Error loading annotated file {atac_annotated_file}: {e}")
                continue

            # Ensure required columns exist in annotated file
            required_cols_annotated = {"ensembl_gene_id"}
            if not required_cols_annotated.issubset(atac_annotated_df.columns):
                print(f"Skipping {atac_annotated_file} - Missing required columns in annotated file.")
                continue

            # Add atac_intersect column to annotated file
            atac_annotated_df.insert(0, "atac_intersect", 0)  # Initialize with 0

            # Fill in "1" for rows where ensembl_gene_id matches gene_id from intersected regions
            for idx, row in atac_annotated_df.iterrows():
                for (chrom, start, end), gene_id in intersected_regions.items():
                    if row["ensembl_gene_id"] == gene_id:
                        atac_annotated_df.at[idx, "atac_intersect"] = 1
                        break

            # Save updated annotated file
            updated_filename = os.path.join(atac_annotated_path, f"updated_{os.path.basename(atac_annotated_file)}")
            atac_annotated_df.to_csv(updated_filename, index=False)
            print(f"Saved updated file: {updated_filename}")




