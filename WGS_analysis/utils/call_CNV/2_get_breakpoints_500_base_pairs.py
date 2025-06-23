import pandas as pd
import numpy as np
import os

# path to project directory
BASE_DIR = "/data6/johnathan/cnv_calling_input"

# CNV calls with different bin sizes obtained from CNVnator
# BIN_SIZES = [50000, 25000, 10000, 5000, 2500, 1000, 500, 200, 100]

# Path to the directory containing CNV calls for each bin size
# CNV_DIR = {bin_size: f"{BASE_DIR}/2_call_cnvs/data/calls_{bin_size}" for bin_size in BIN_SIZES}

# File containing list of samples
SAMPLES_FILE = f"{BASE_DIR}/slurm/files/0_smap.txt"

# details for chromosome of interest
CHR = "chr18"
CHR_START = 1
CHR_END = 80373300
EXCLUDE_START = 200000
EXCLUDE_END = 200000
EXCLUDE_CENTROMERE_START = 15000000
EXCLUDE_CENTROMERE_END = 21000000

# CNV type of interest
CNV_TYPE = "duplication"

# Path to the CSV file to save breakpoints
BREAKPOINT_CSV = f"/path/wgs_processing_{CNV_TYPE}_500bp_breakpoints_with_gene_name_auto.csv"

# function to find nearest value in an array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

gtf_file = "/path/gene_annotations/gencode.v44.basic.annotation.gtf"

gene_df = pd.read_csv(gtf_file, usecols=[0,2,3,4,8], names=["chrm", "product_type", "start", "end", "product_info"], header=None, sep="\t", comment="#")
gene_df = gene_df.loc[gene_df.product_type=="gene"]
gene_df["gene_length"] = gene_df.end-gene_df.start
# create gene id, type and name columns
gene_df["gene_id"] = gene_df.product_info.apply(lambda x: x.split(";")[0].strip("gene_id ").strip('"').split(".")[0])
gene_df["gene_type"] = gene_df.product_info.apply(lambda x: x.split(";")[1].strip("gene_type ").strip('"'))
gene_df["gene_name"] = gene_df.product_info.apply(lambda x: x.split(";")[2].strip("gene_name ").strip('"'))
# gene_df = gene_df[gene_df['chrm'] == 'chr18']

# get sample IDs
# f = open(SAMPLES_FILE, "r")
# samples = f.read().splitlines()
# f.close()

samples = open('/path/slurm/files/0_smap.txt', 'r').read().splitlines()

breakpoints = []
for sample in samples:
    print(sample)
    # Start with the largest bin size and go down to 
    # the smallest bin size until a CNV is found
    # read in CNV calls
    cnv_file = f"/path//calls/{sample}_500.tsv"
    cnv_df = pd.read_csv(cnv_file, sep="\t")

    df = cnv_df[cnv_df["CNV_type"] == CNV_TYPE]

    for index, row in df.iterrows():
        matching_lines = gene_df[(gene_df['chrm'] == row["chrm"]) & (gene_df['gene_type']=='protein_coding') & (((gene_df.start>=row.start) & (gene_df.start<=row.end)) | ((gene_df.end>=row.start) & (gene_df.end<=row.end)))]
        gene_list = '|'.join(matching_lines.gene_name.to_list())
        gene_id_list = '|'.join(matching_lines.gene_id.tolist())
        
        # save the start and end positions of all CNV calls
        breakpoint = {
            "sample": sample,
            "chr": row["chrm"],
            "start": row["start"],
            "end": row["end"],
            "bin_size": 500,
            "homozygous": row["read_depth"] < 0.25,
            "gene_name": gene_list,
            "gene_id": gene_id_list
        }
        breakpoints.append(breakpoint)

breakpoints_df = pd.DataFrame(breakpoints)

breakpoints_df.to_csv(BREAKPOINT_CSV, index=False)
