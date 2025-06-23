import pandas as pd
import scipy.stats as stats
from scipy import stats

# Change these variables for different runs
deg_filename="/Users/jks6575/Dropbox/Jiawan/Analysis/outlier_expression/DEGs/CRISPR_DEGs/mn_CRISPR_DEG.csv"
output_filename="/Users/jks6575/Dropbox/Jiawan/Analysis/outlier_expression/DEGs/mn_CRISPR_fishers_exact_coding_noncoding.csv"
samples=['CR001_A38', 'CR001_WT']

# Load in variant data
vardf=pd.read_csv('/Users/jks6575/Dropbox/Jiawan/Analysis/outlier_expression/variant_data/ipsc_family_variants_updated.csv')
vardf = vardf[vardf['Sample'].isin(samples)]
# filter variant in CR001_A38 the same with CR001_WT
vardf = vardf[
    (vardf['Sample'] == 'CR001_A38') & 
    (vardf['vid'].isin(vardf[vardf['Sample'] == 'CR001_WT']['vid']))
]
print(vardf)
vardf.to_csv("A38_overlap_variant.csv", index=False)
variant_types=[ 'snv_coding','snv_noncoding', 'deletion_coding', 'deletion_noncoding', 'duplication_coding', 
               'duplication_noncoding', 'STR_coding', 'STR_noncoding'] #'SNV_coding',
# Load in gene list
degs=pd.read_csv(deg_filename).ensembl_gene_id.to_list()
print(degs)

# Perfrom Fisher's Exact tests separately on over and under expressed genes in each cell type on specific sets of samples
# For example, consider samples with both RNA-seq and WGS data
stat_lst=[]
gene_df = pd.read_csv('/Users/jks6575/Dropbox/Jiawan/Analysis/outlier_expression/DEGs/mature_neuron_counts_tpm.csv')
gene_df['Gene'] = gene_df['Gene'].astype(str)
gene_space = sorted(list(gene_df['Gene'].unique()))

for vt in variant_types:
    df=pd.DataFrame({'Gene':gene_space})
    print(df)
    # Check for variant
    df['Variant']=0
    df.loc[df.Gene.isin(vardf[(vardf.variant_type_2==vt)].gene_id.to_list()), 'Variant']=1
    print(df)
    
    # Check for DEG
    df['DEG']=0
    df.loc[df.Gene.isin(degs), 'DEG']=1
    
    # Fisher's exact test for if a gene (1) has outlier expression in a sample and (2) has a variant in the same sample
    count_df=df[['Variant', 'DEG']].groupby(['Variant', 'DEG']).size().to_frame()
    print(count_df)
    count_df.reset_index(inplace=True)
    count_df=count_df.pivot(index='DEG', columns='Variant', values=0)
    count_df.fillna(0, inplace=True)
    count_df=count_df.astype(int)
    res=stats.fisher_exact(count_df)
    or_res=stats.contingency.odds_ratio(count_df)
    ci=or_res.confidence_interval()
    stat_lst.append([vt, res.statistic, ci[0], ci[1], res.pvalue])
statdf=pd.DataFrame(stat_lst, columns=['Variant', 'Odds ratio', '95% C.I. lower', '95% C.I. upper', 'p value'])

# Save to file
statdf['BH FDR']=stats.false_discovery_control(statdf['p value'].to_numpy(), method='bh')

print(statdf[statdf['BH FDR']<=0.05])

# Save to file
statdf.to_csv(output_filename, index=False)
