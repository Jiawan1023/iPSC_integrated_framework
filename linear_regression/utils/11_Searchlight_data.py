import pandas as pd

import os
cwd = os.getcwd()
dropbox=cwd.split('Dropbox')[0]+'Dropbox/'
dropbox=dropbox.replace('\\', '/')

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Use data from WGS paper
data_loc=dropbox+"16p12.2 project/Human patients project/WGS paper/CLEAN_CODE/Searchlight/Searchlight_Data"
svip=pd.DataFrame()
for cohort in ['deletion', 'duplication']:
	cdf=pd.read_csv(f'{data_loc}/16p11.2 {cohort}.csv')
	svip=pd.concat([svip, cdf])

# Add in additional data
summ_var=pd.read_csv(f'{dropbox}/SVIP Project/Phenotype information/SVIP_2019_phenotypic_data/Simons_Searchlight_Phase1_16p11.2_Dataset_v11.0/svip_summary_variables.csv')
summ_var['Sample']=summ_var.sfari_id.str.replace('-', '.')
summ_var.index=summ_var.Sample.to_list()
svip['Age']=svip.Sample.map(summ_var.age_months.to_dict())

# Restrict to just the needed data
svip['SNV']=True
svip.loc[svip['Missense'].isnull(), 'SNV']=False
svip['CNV']=True
svip.loc[svip['Genes del.'].isnull(), 'CNV']=False
svip=svip[['Sample', 'Cohort', 'SNV', 'CNV', 'Sex', 'Age',
			'Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Autism behavior (BSI)', 'BMI z-score', 'Head circumference z-score']]
svip.sort_values(by='Sample', inplace=True)

# Save
svip.to_csv('files/Searchlight_phenotypes.csv', index=False)

# Make gene lists
wgs=dropbox+"16p12.2 project/Human patients project/WGS paper/"

# SNVs
snv=pd.read_csv(wgs+'32_SVIP analysis/SNVs_indels/SVIP_rare_deleterious_snvs_indels.csv')
snv=snv[snv.Sample.isin(svip.Sample.to_list())]
snv['Gene_id_']=snv['Gene_id_'].str.split(';')
snv['Gene_symbol']=snv.Gene_symbol.str.split(';')
snv=snv.explode(['Gene_id_', 'Gene_symbol'])
snv=snv[['Sample', 'Gene_id_', 'Gene_symbol']]
snv['Variant']='SNV'

# CNVs
cnv=pd.read_csv(wgs+'32_SVIP analysis/CNVs/SVIP_CNVs_by_gene_loeuf.csv')
cnv=cnv[cnv.Sample.isin(svip.Sample.to_list())]
cnv=cnv[cnv.CNV_Type.isin(['Del', 'Dup'])]
cnv['Gene_id_']=cnv['Gene_id_'].str.split(';')
cnv['Gene']=cnv.Gene.str.split(';')
cnv=cnv.explode(['Gene_id_', 'Gene'])
cnv=cnv[['Sample', 'Gene_id_', 'Gene']]
cnv.columns=['Sample', 'Gene_id_', 'Gene_symbol']
cnv['Variant']='CNV'

# Save
svip_genes=pd.concat([snv, cnv])
svip_genes.sort_values(by=['Sample', 'Variant', 'Gene_symbol'], inplace=True)
svip_genes.to_csv('files/Searchlight_genes.csv', index=False)

# Load in module and pathway gene sets
# For each proband, count the number of genes they have in the set
genesets=[i for i in os.listdir(dropbox+'Jiawan/Analysis/WGCNA/Brainspan/module_genes') if '.csv' in i]

svip_samps=sorted(list(svip_genes.Sample.unique()))
outdf=pd.DataFrame(index=svip_samps)
for gs in genesets:
	print(gs)
	gsdf=pd.read_csv(dropbox+"Jiawan/Analysis/WGCNA/Brainspan/module_genes/"+gs)
	name=gs.split('.csv')[0]
	if 'genenames' in name:
		name=name.split('_')[0]
	
	if 'gene' in gsdf.columns.to_list():
		gene_list=gsdf.gene.to_list()
	else:
		gene_list=gsdf.ensembl_gene_id.to_list()
	
	gs_svip=svip_genes[svip_genes.Gene_id_.isin(gene_list)]
	gs_samps=gs_svip.Sample.value_counts()
	count_df=pd.DataFrame(gs_samps)
	count_df.columns=[name]
	outdf=pd.merge(outdf, count_df, right_index=True, left_index=True, how='outer')
	
outdf.fillna(0, inplace=True)
outdf=outdf.astype(int)

pdf=PdfPages('Figures/11_Searchlight_gene_set_counts.pdf')
for col in outdf.columns.to_list():
	sns.histplot(data=outdf, x=col, discrete=True)
	plt.title(col)
	pdf.savefig()
	plt.close()
pdf.close()

# Save
outdf['Sample']=outdf.index.to_list()
outdf.to_csv('files/Searchlight_gene_set_counts.csv', index=False)
