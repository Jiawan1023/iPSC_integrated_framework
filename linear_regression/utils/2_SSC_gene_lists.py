import pandas as pd

import os
cwd = os.getcwd()
dropbox=cwd.split('Dropbox')[0]+'Dropbox/'
dropbox=dropbox.replace('\\', '/')
wgs=dropbox+"16p12.2 project/Human patients project/WGS paper/"

# Load in samples
ssc_samps=pd.read_csv('files/SSC_phenotypes.csv').Sample.to_list()

# Create lists of affected genes for each sample
# SNVs
snv=pd.read_csv(wgs+'31_SSC analysis/SNV_indel_variants/SSC_rare_deleterious_snvs_indels_dbd.csv')
snv=snv[snv.Sample.str.contains('.p1')]
snv=snv[snv.Sample.isin(ssc_samps)]
# Remove SNV first hits
dbd_fh=pd.read_csv(wgs+"31_SSC analysis/Genotype-phenotype integration/1_variant_preparation/first_hit_variants/dbd_tier1_snvs.csv")
dbd_fh['Sample_Gene']=dbd_fh.Sample+'.'+dbd_fh.Gene_id_
snv['Sample_Gene']=snv.Sample+'.'+snv.Gene_id_
# If sample has multiple variants in the same first hit gene, consider the one with the highest CADD score as the first hit
# If they have the same CADD, take the first one
snv['First_hit']=False
snv.loc[snv.Sample_Gene.isin(dbd_fh.Sample_Gene.to_list()), 'First_hit']=True
repeats=snv[snv.Sample_Gene.isin(dbd_fh.Sample_Gene.to_list())]['Sample_Gene'].value_counts()
repeats=repeats[repeats>1].index.to_list()
for rp in repeats:
	snv.loc[snv.Sample_Gene==rp, 'First_hit']=False
	cadds=[i for i in snv[snv.Sample_Gene==rp]['CADD_PHRED'].to_list() if i==i]
	if len(cadds)>0:
		idx=snv[(snv.Sample_Gene==rp) & (snv.CADD_PHRED==max(cadds))].index.to_list()
	else:
		idx=snv[(snv.Sample_Gene==rp)].index.to_list()
	snv.loc[idx[0], 'First_hit']=True
snv=snv[~snv.First_hit]
snv['Gene_id_']=snv['Gene_id_'].str.split(';')
snv['Gene_symbol']=snv.Gene_symbol.str.split(';')
snv=snv.explode(['Gene_id_', 'Gene_symbol'])
snv=snv[['Sample', 'Gene_id_', 'Gene_symbol']]
snv['Variant']='SNV'

# CNVs
cnv=pd.read_csv(wgs+'31_SSC analysis/CNVs/SSC_Sanders_CNVs_by_gene_loeuf.csv')
cnv=cnv[cnv.Sample.str.contains('.p1')]
cnv=cnv[cnv.Sample.isin(ssc_samps)]
cnv=cnv[cnv.CNV_Type.isin(['Del', 'Dup'])]
# Remove first hit CNVs
for cohort in ['large_rare_deletions', 'large_rare_duplications']:
	fh=pd.read_csv(wgs+"31_SSC analysis/Genotype-phenotype integration/1_variant_preparation/first_hit_variants/"+cohort+'.csv')
	fh['Gene_id_']=fh.Gene_id_.str.split('.', expand=True)[0]
	cnv_type='Del'
	if cohort=='large_rare_duplications':
		cnv_type='Dup'
	fh['Sample_gene_type']=fh.Sample+'.'+fh.Gene_id_+'.'+cnv_type

	cnv['Sample_gene_type']=cnv.Sample+'.'+cnv.Gene_id_+'.'+cnv.CNV_Type

	cnv=cnv[~cnv.Sample_gene_type.isin(fh['Sample_gene_type'].to_list())]

cnv['Gene_id_']=cnv['Gene_id_'].str.split(';')
cnv['Gene']=cnv.Gene.str.split(';')
cnv=cnv.explode(['Gene_id_', 'Gene'])
cnv=cnv[['Sample', 'Gene_id_', 'Gene']]
cnv.columns=['Sample', 'Gene_id_', 'Gene_symbol']
cnv['Variant']='CNV'

# STRs
str=pd.read_csv(wgs+'31_SSC analysis/STR_variants/SSC_Exonic_STRs.csv')
str=str[str.Sample.str.contains('.p1')]
str=str[str.Sample.isin(ssc_samps)]

str['Gene_id_']=str['Gene_id_'].str.split(';')
str['Gene_symbol']=str.Gene_symbol.str.split(';')
str=str.explode(['Gene_id_', 'Gene_symbol'])
str=str[['Sample', 'Gene_id_', 'Gene_symbol']]
str['Variant']='STR'

df=pd.concat([snv, cnv, str])
df.sort_values(by=['Sample', 'Variant', 'Gene_symbol'], inplace=True)
df.to_csv('files/SSC_genes.csv', index=False)
