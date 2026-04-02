import pandas as pd

# Reformat CNV calls

# Input and output
CNV='../data/filtered_calls/iPSC_cohort_calls.csv'
OUTPUT='../data/filtered_calls/CNV_calls.csv'

annos=['upstream', 'utr5', 'intron', 'exon', 'utr3', 'downstream']

# Load data
df=pd.read_csv(CNV)

samp_cols=[i for i in df.columns.to_list() if 'SG' in i or 'CR' in i]

df=df[['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq']+samp_cols+annos]

# Pivot
df=df.melt(id_vars=['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq']+samp_cols,
		   	var_name='Var_type', value_name='gene_id')
df=df[~df.gene_id.isnull()]

# Explode by Gene ID
df.gene_id=df.gene_id.str.split(';')
df=df.explode('gene_id')

# Rename variant type
df['variant_type']=df.CNV_type+'_'+df.Var_type
df['c_nc']='noncoding'
df.loc[df.Var_type=='exon', 'c_nc']='coding'
df['variant_type_2']=df.CNV_type+'_'+df.c_nc

df['vid']=df.Chr+'_'+df.Start.astype(str)+'_'+df.End.astype(str)+'_'+df.CNV_type

# Save
df=df[['Sample', 'gene_id', 'vid', 'variant_type', 'variant_type_2']+samp_cols]
df.to_csv(OUTPUT, index=False)

print(df)
