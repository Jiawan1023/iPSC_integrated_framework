import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype']=42

# Plot regression results
outdf=pd.read_csv('Result_tables/9_16p12_regression.csv')
outdf['Phenotype']=outdf.Phenotype.str.replace(' (Child domain)', '')

wgcna=['blue', 'brown', 'green', 'magenta', 'turquoise']
pathway=['BMP', 'MAPK_ERK', 'NOTCH', 'RAS', 'RHO', 'WNT']

quant_phenos=['De Vries Score', 'BMI Z Score', 'Head Circumference Z Score', 'HRS-MAT', 'SRS Raw Score']

# Plot results
pdf=PdfPages('Figures/10_16p12_results.pdf')
for phenos in [quant_phenos]:
	sup_idx=[quant_phenos].index(phenos)
	sup_title=['Binary phenotypes', 'Quantitative phenotypes'][sup_idx]
	fig, axs=plt.subplots(ncols=2, figsize=(10, 5))
	for gs in [wgcna, pathway]:
		idx=[wgcna, pathway].index(gs)
		
		plotdf=outdf[(outdf.Phenotype.isin(phenos)) & (outdf.Variable.isin(gs))].copy()
		plotdf['star']=''
		plotdf.loc[plotdf['p value']<=0.05, 'star']='*'
		plotdf.loc[plotdf['BH FDR']<=0.05, 'star']='**'
		
		heatdf=plotdf.pivot(index='Gene set', columns='Phenotype', values='Estimate')
		stardf=plotdf.pivot(index='Gene set', columns='Phenotype', values='star')
		
		heatdf=heatdf.loc[gs, phenos]
		stardf=stardf.loc[gs, phenos]
		
		biggest=max(abs(plotdf["Estimate"].to_numpy()))
		
		colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
		cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
		cmap.set_bad('#CCCCCC')
		sns.heatmap(data=heatdf, cmap=cmap, vmin=-biggest, vmax=biggest, square=True, fmt='', linecolor='k', linewidths=0.75, annot=stardf, ax=axs[idx])
		
		title=['WGCNA', 'Pathway'][idx]
		axs[idx].set_title(title)
	
	plt.suptitle(sup_title)
	plt.tight_layout()
	pdf.savefig()
	plt.close()
pdf.close()
		
