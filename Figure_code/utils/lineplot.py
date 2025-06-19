import pandas as pd

# For plotting
import seaborn as sns
import matplotlib.pyplot as plt

# Get just the genes you want
# Make a dictionary of gene sets (list of genes + name) examples
gene_sets={'pluripotency_marker':['ENSG00000181449', 'ENSG00000204531', 'ENSG00000196371'],
			'npc_marker':['ENSG00000181449', 'ENSG00000114315', 'ENSG00000132688', 'ENSG00000170370'],
			'neuron_marker':['ENSG00000078018', 'ENSG00000132535', 'ENSG00000149294', 'ENSG00000167281']}

# Load in the normalized data as df
# Transpose the data (its easier to work with if samples are rows and not columns)
tdf=df.transpose()
# For each gene set defined above, get the average normalized expression for each row
for gs in gene_sets.keys():
	# Check that all the genes you want are in the data
	gl=gene_sets[gs]
	gl_in=list(set(gl).intersection(set(tdf.columns.to_list())))
	gl_out=list(set(gl)-set(gl_in))
	if len(gl_out)>0:
		print('Genes in', gs, 'not found in data:', ','.join(gl_out))
	# Calculate average for genes found in data
	tdf[gs]=tdf[gl_in].mean(axis=1)
# Clean up the dataframe by only keeping the values we want for plotting
tdf=tdf[list(gene_sets.keys())]

# Convert data to long form
tdf['Sample']=tdf.index.to_list()
ldf=pd.melt(tdf, id_vars='Sample', var_name='gene_set', value_name='average_expression')

# Annotate rows with cell type
ldf['celltype']=''
ldf.loc[ldf.Sample.str.lower().str.contains('ipsc'), 'celltype']='iPSC'
ldf.loc[ldf.Sample.str.lower().str.contains('npc'), 'celltype']='NPC'
ldf.loc[ldf.Sample.str.contains('IM'), 'celltype']='IM'
ldf.loc[ldf.Sample.str.contains('MN'), 'celltype']='MN'
ldf.celltype=pd.Categorical(ldf.celltype, ['iPSC', 'NPC', 'IM', 'MN'])
ldf.sort_values(by='celltype', inplace=True)
ldf

# Plot gene set expression over celltypes
sns.lineplot(data=ldf, x='celltype', y='average_expression', hue='gene_set', palette="Dark2")
plt.legend(bbox_to_anchor=(1.45, 0.5), loc='right')
my_path= "#your_path"
plt.savefig(my_path + 'Marker_lineplot.pdf')
