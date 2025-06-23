import pandas as pd

import os
cwd = os.getcwd()
dropbox=cwd.split('Dropbox')[0]+'Dropbox/'
dropbox=dropbox.replace('\\', '/')

# Load in module and pathway gene sets
# For each proband, count the number of genes they have in the set
genesets=['input_gene_set_file_names']

ssc_genes=pd.read_csv('files/SSC_genes.csv')
ssc_samps=sorted(list(ssc_genes.Sample.unique()))
outdf=pd.DataFrame(index=ssc_samps)
for gs in genesets:
	gsdf=pd.read_csv(dropbox+f"Jiawan/Analysis/big_cohort_application/{gs}.csv")
	name=gs.split('_g')[0]
	
	if 'gene' in gsdf.columns.to_list():
		gene_list=gsdf.gene.to_list()
	else:
		gene_list=gsdf.ensembl_gene_id.to_list()
	
	gs_ssc=ssc_genes[ssc_genes.Gene_id_.isin(gene_list)]
	gs_samps=gs_ssc.Sample.value_counts()
	count_df=pd.DataFrame(gs_samps)
	count_df.columns=[name]
	outdf=pd.merge(outdf, count_df, right_index=True, left_index=True, how='outer')
	
outdf.fillna(0, inplace=True)
outdf=outdf.astype(int)

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf=PdfPages('Figures/3_SSC_gene_set_counts.pdf')
for col in outdf.columns.to_list():
	sns.histplot(data=outdf, x=col)
	plt.title(col)
	pdf.savefig()
	plt.close()
pdf.close()

# Save
outdf['Sample']=outdf.index.to_list()
outdf.to_csv('files/SSC_gene_set_counts.csv', index=False)
