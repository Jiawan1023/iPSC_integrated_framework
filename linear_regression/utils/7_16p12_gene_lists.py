import pandas as pd
import numpy as np

import os
cwd = os.getcwd()
dropbox=cwd.split('Dropbox')[0]+'Dropbox/'
dropbox=dropbox.replace('\\', '/')
wgs=dropbox+"16p12.2 project/Human patients project/WGS paper/"

# Load in samples
samps=pd.read_csv('files/16p12_phenotypes.csv').Sample.to_list()

# SNVs
snvs=pd.read_csv(wgs+"9_Rare variant calls/Rare coding SNV calls/Rare_Deleterious_Exonic_Variants.csv")
snvs=snvs[snvs.Sample.isin(samps)]
snvs['Gene_id_']=snvs['Gene_id_'].str.split(';')
snvs=snvs.explode('Gene_id_')
snvs['Variant']='SNV'

# CNVs
cnvs=pd.read_csv(wgs+"7_Structural variants/sv_calls_combined.txt", sep='\t')
cnvs=cnvs[cnvs.Sample.isin(samps)]
cnvs['Gene_ID']=cnvs['Gene_ID'].str.split(' ')
cnvs=cnvs.explode('Gene_ID')
cnvs.replace('.', np.nan, inplace=True)
cnvs=cnvs[['Sample', 'Gene_ID', 'Gene_Symbol']]
cnvs.columns=['Sample', 'Gene_id_', 'Gene_symbol']
cnvs['Variant']='CNV'

# STRs
strs=pd.read_csv(wgs+"8_STR variants/Exonic_STR_expansions.csv", index_col=0)
strs=strs[strs.Sample.isin(samps)]
strs['Gene_id_']=strs['Gene_id_'].str.split(';')
strs=strs.explode('Gene_id_')
strs['Variant']='STR'

df=pd.concat([snvs, cnvs, strs])
df.sort_values(by=['Sample', 'Variant', 'Gene_symbol'], inplace=True)
df=df[['Sample', 'Gene_id_', 'Gene_symbol', 'Variant']]
df.to_csv('files/16p12_genes.csv', index=False)
