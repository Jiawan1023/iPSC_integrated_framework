import pandas as pd

import os
cwd = os.getcwd()
dropbox=cwd.split('Dropbox')[0]+'Dropbox/'
dropbox=dropbox.replace('\\', '/')

# Load in module and pathway gene sets
# For each proband, count the number of genes they have in the set
genesets=['input_module_pathway_gene_name']

genes=pd.read_csv('files/16p12_genes.csv')
samps=sorted(list(genes.Sample.unique()))
outdf=pd.DataFrame(index=samps)
for gs in genesets:
	gsdf=pd.read_csv(dropbox+f"Jiawan/Analysis/big_cohort_application/{gs}.csv")
	name=gs.split('_g')[0]
	
	if 'gene' in gsdf.columns.to_list():
		gene_list=gsdf.gene.to_list()
	else:
		gene_list=gsdf.ensembl_gene_id.to_list()
	
	gs_16p12=genes[genes.Gene_id_.isin(gene_list)]
	gs_samps=gs_16p12.Sample.value_counts()
	count_df=pd.DataFrame(gs_samps)
	count_df.columns=[name]
	outdf=pd.merge(outdf, count_df, right_index=True, left_index=True, how='outer')
	
outdf.fillna(0, inplace=True)
outdf=outdf.astype(int)

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf=PdfPages('Figures/8_16p12_gene_set_counts.pdf')
for col in outdf.columns.to_list():
	sns.histplot(data=outdf, x=col)
	plt.title(col)
	pdf.savefig()
	plt.close()
pdf.close()

# Save
outdf['Sample']=outdf.index.to_list()
outdf.to_csv('files/16p12_gene_set_counts.csv', index=False)
