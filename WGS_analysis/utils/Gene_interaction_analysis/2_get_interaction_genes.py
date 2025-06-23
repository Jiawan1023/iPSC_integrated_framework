import pandas as pd
import numpy as np
import os
from scipy.stats import f_oneway
from statsmodels.formula.api import ols
import statsmodels.api as sm

input_dir = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/imn/gene/fam3"
output_dir = "/Users/jks6575/Dropbox/Jiawan/Analysis/Genetic_interaction_new/imn/gene"


# Load TPM counts file (assumed to be the same for all variants)
ipsc_counts_tpm = pd.read_csv(os.path.join(input_dir, "immature_counts_tpm.csv"))
# Extract relevant columns
columns_to_extract = [
    "Gene",
    "X3110.IM.R1", "X3110.IM.R2", "X3110.IM.R3", "X321.IM.R1", "X321.IM.R2", "X321.IM.R3",
    "X331.IM.R1", "X331.IM.R2", "X331.IM.R3", "X342.IM.R1", "X342.IM.R2", "X342.IM.R3"
]
ipsc_counts_tpm = ipsc_counts_tpm[columns_to_extract]


def parse_groups_from_filename(filename):

    # Remove file extension and split by underscores
    parts = filename.replace("_imn.csv", "").split("_")

    # Extract groups using "versus" as separator
    try:
        # Identify indices where "versus" appears
        versus_indices = [i for i, part in enumerate(parts) if part == "versus"]

        # Extract groups based on positions in filename
        group_A = parts[1:versus_indices[0]]  # Numbers before first 'versus'
        group_B = parts[versus_indices[0] + 1 : versus_indices[1]]  # Between first and second 'versus'
        group_C = parts[versus_indices[1] + 1 : versus_indices[2]]  # Between second and third 'versus'
        group_D = parts[versus_indices[2] + 1 :]  # After third 'versus'
    
    except IndexError:
        raise ValueError(f"Filename format is incorrect: {filename}")

    # Convert to expected file names
    group_A = [f"X{v}.IM.R{i}" for v in group_A for i in range(1, 4)]
    group_B = [f"X{v}.IM.R{i}" for v in group_B for i in range(1, 4)]
    group_C = [f"X{v}.IM.R{i}" for v in group_C for i in range(1, 4)]
    group_D = [f"X{v}.IM.R{i}" for v in group_D for i in range(1, 4)]

    return group_A, group_B, group_C, group_D


# Function to run ANOVA
def run_anova(row, group_A, group_B, group_C, group_D):
    """
    Runs a two-way ANOVA test on the given row, ensuring that only existing columns are used.
    """
    # Combine all groups
    all_groups = group_A + group_B + group_C + group_D

    # Filter out columns that are not in the row index (i.e., skip missing columns)
    existing_groups = [col for col in all_groups if col in row.index]

    if not existing_groups:
        return {"deletion_p_value": None, "variant_p_value": None, "interaction_p_value": None}  # Skip if no valid data

    data = {
        "group1": (["deletion"] * len(group_A) + 
                   ["nondeletion"] * len(group_B) + 
                   ["nondeletion"] * len(group_C) + 
                   ["deletion"] * len(group_D))[:len(existing_groups)],  # Trim to match existing columns

        "value": row[existing_groups].values.flatten(),  # Use only existing columns

        "group2": (["variant"] * len(group_A) + 
                   ["variant"] * len(group_B) + 
                   ["nonvariant"] * len(group_C) + 
                   ["nonvariant"] * len(group_D))[:len(existing_groups)],  # Trim to match existing columns
    }

    row_df = pd.DataFrame(data)
    row_df["value"] = pd.to_numeric(row_df["value"], errors="coerce")

    # Fit ANOVA model
    model = ols('value ~ group1 * group2', data=row_df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    return {
        "deletion_p_value": anova_table.loc["group1", "PR(>F)"] if "group1" in anova_table.index else None,
        "variant_p_value": anova_table.loc["group2", "PR(>F)"] if "group2" in anova_table.index else None,
        "interaction_p_value": anova_table.loc["group1:group2", "PR(>F)"] if "group1:group2" in anova_table.index else None,
    }


# Process all variant files
for file in os.listdir(input_dir):
    if file.startswith("variant_") and file.endswith(".csv"):
        print(f"Processing {file}...")
        
        # Load variant file
        variant_df = pd.read_csv(os.path.join(input_dir, file))
        filtered_data = ipsc_counts_tpm[ipsc_counts_tpm['Gene'].isin(variant_df['gene_id'])]

        # Extract dynamic groups
        group_A, group_B, group_C, group_D = parse_groups_from_filename(file)

        # Run ANOVA
        anova_results = filtered_data.apply(lambda row: run_anova(row, group_A, group_B, group_C, group_D), axis=1, result_type="expand")

        # Merge results
        final_results = pd.concat([filtered_data, anova_results], axis=1)
        significant_results = final_results[final_results["interaction_p_value"] <= 0.05]

        # Save output
        output_file = os.path.join(output_dir, f"significant_results_{file}")
        significant_results.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
