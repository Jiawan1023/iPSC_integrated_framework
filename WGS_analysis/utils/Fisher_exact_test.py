import pandas as pd
import scipy.stats as stats
from scipy import stats

# Change these variables for different runs

output_filename="Fisher_exact_test/remove/mn_CRISPR_fishers_exact_coding_noncoding_remove.csv"

vardf=pd.read_csv('isogenic_overlap_variant.csv')
unique_df = pd.read_csv('isogenic_unique_variant.csv')
deg_df = pd.read_csv('Fisher_exact_test/remove/common_effect_gene_expression/mn_isogenic_DEGs_with_variant.csv')

#create DEG table, with annotation of overlap between genes and variants
# Make gene sets for fast lookup
shared_gene_set = set(vardf['gene_id'].dropna())
unique_gene_set = set(unique_df['gene_id'].dropna())

# Create new columns
deg_df['shared_rare_variant'] = deg_df['ensembl_gene_id'].isin(shared_gene_set).astype(int)
deg_df['unique_rare_variant'] = deg_df['ensembl_gene_id'].isin(unique_gene_set).astype(int)

# (Optional) save updated dataframe
deg_df.to_csv(
    'Fisher_exact_test/remove/mn_isogenic_DEGs_with_variant_updated.csv',
    index=False
)

import numpy as np

#remove genes from DEG overlap with unique variants or common variants
# --- 1) Remove rows whose gene is in unique_gene_set (i.e., unique_rare_variant == 1) ---
deg_df = deg_df[~deg_df['ensembl_gene_id'].isin(unique_gene_set)].copy()


# --- 2) Remove rows that meet BOTH:
#   (a) common_variant == 1
#   (b) GTEX_slope AND PsychENCODE_slope are the same direction as log2FoldChange
# ---

def same_direction(a, b):
    """
    True if a and b are both non-missing, non-zero, and have the same sign.
    """
    a = pd.to_numeric(a, errors="coerce")
    b = pd.to_numeric(b, errors="coerce")
    return (np.sign(a) == np.sign(b)) & (a != 0) & (b != 0)

common_is_1 = pd.to_numeric(deg_df['common_variant'], errors="coerce").fillna(0).astype(int).eq(1)

gtex_same = same_direction(deg_df['GTEX_slope'], deg_df['log2FoldChange'])
psy_same  = same_direction(deg_df['PsychENCODE_slope'], deg_df['log2FoldChange'])

remove_mask = common_is_1 & gtex_same & psy_same

deg_df = deg_df[~remove_mask].copy()
variant_types=[ 'snv_coding', 'snv_noncoding', 'deletion_coding', 'deletion_noncoding',  'STR_coding', 'STR_noncoding'] #change to different types,
# Load in gene list
degs=deg_df.ensembl_gene_id.to_list()
print(len(degs))

# Perfrom Fisher's Exact tests separately on over and under expressed genes in each cell type on specific sets of samples
# For example, consider samples with both RNA-seq and WGS data
stat_lst=[]
gene_df = pd.read_csv('C:/Users/Jiawan/Dropbox/Jiawan/Analysis/outlier_expression/DEGs/mn_counts_tpm.csv')
gene_df['Gene'] = gene_df['Gene'].astype(str)
gene_space = sorted(list(gene_df['Gene'].unique()))

for vt in variant_types:
    df=pd.DataFrame({'Gene':gene_space})
    # Check for variant
    df['Variant']=0
    df.loc[df.Gene.isin(vardf[(vardf.variant_type_2==vt)].gene_id.to_list()), 'Variant']=1
    
    # Check for DEG
    df['DEG']=0
    df.loc[df.Gene.isin(degs), 'DEG']=1
    
    
    # Fisher's exact test for if a gene (1) has outlier expression in a sample and (2) has a variant in the same sample
    count_df=df[['Variant', 'DEG']].groupby(['Variant', 'DEG']).size().to_frame()
    count_df.reset_index(inplace=True)
    count_df=count_df.pivot(index='DEG', columns='Variant', values=0)
    count_df.fillna(0, inplace=True)
    count_df=count_df.astype(int)
    print(count_df)
    res=stats.fisher_exact(count_df)
    or_res=stats.contingency.odds_ratio(count_df)
    ci=or_res.confidence_interval()
    stat_lst.append([vt, res.statistic, ci[0], ci[1], res.pvalue])
statdf=pd.DataFrame(stat_lst, columns=['Variant', 'Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value'])

# Save to file
#statdf['BH FDR']=stats.false_discovery_control(statdf['p value'].to_numpy(), method='bh')

#print(statdf[statdf['BH FDR']<=0.05])
from statsmodels.stats.multitest import multipletests

statdf['BH FDR'] = multipletests(
    statdf['p value'].to_numpy(),
    method='fdr_bh'
)[1]

print(statdf[statdf['BH FDR'] <= 0.05])

# Save to file
statdf.to_csv(output_filename, index=False)

import pandas as pd
from scipy import stats

stat_lst = []

df = pd.DataFrame({"Gene": gene_space})

# Check for variant
df["Variant"] = 0
df.loc[df["Gene"].isin(vardf["gene_id"].astype(str)), "Variant"] = 1

# Check for DEG
df["DEG"] = 0
df.loc[df["Gene"].isin(pd.Series(degs).astype(str)), "DEG"] = 1

# Build robust 2x2 table: rows=DEG(0/1), cols=Variant(0/1)
count_df = (
    df.groupby(["DEG", "Variant"])
      .size()
      .unstack("Variant", fill_value=0)
      .reindex(index=[0, 1], columns=[0, 1], fill_value=0)
      .astype(int)
)

print(count_df)

table = count_df.to_numpy()   # guaranteed shape (2,2)

res = stats.fisher_exact(table)

or_res = stats.contingency.odds_ratio(table)
ci = or_res.confidence_interval()

stat_lst.append([res.statistic, ci.low, ci.high, res.pvalue])

print(res)

statdf_all = pd.DataFrame(stat_lst, columns=["Odds ratio", "95% CI low", "95% CI high", "p value"])
print(statdf_all)

output_filename_2="Fisher_exact_test/remove/mn_CRISPR_fishers_exact_all_remove.csv"

statdf_all.to_csv(output_filename_2, index=False)
