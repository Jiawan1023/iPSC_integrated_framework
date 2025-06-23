import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cmx
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype']=42

# Plot regression results

# Load in data
outdf=pd.read_csv('Result_tables/4_SSC_regression.csv')

quant_phenos=['Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score',  'HC_zscore']
primary_variants=['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'Any variant']
wgcna=['input_module_names']
pathway=['input_pathway_names']

# Plot results
for var in ['Interaction', 'Geneset']:
	pdf=PdfPages(f'Figures/5_regression_{var}.pdf')
	for gs in [wgcna, pathway]:
		if var=='Interaction':
			plotdf=outdf[(outdf.Variable=='Interaction') & (outdf['Gene set'].isin(gs))].copy()
		else:
			plotdf=outdf[(outdf.Variable.isin(gs)) & (outdf['Gene set'].isin(gs))].copy()

		plotdf['Gene set']=pd.Categorical(plotdf["Gene set"], gs)
		plotdf['Primary variant']=pd.Categorical(plotdf['Primary variant'], primary_variants)
		plotdf['Phenotype']=pd.Categorical(plotdf['Phenotype'], quant_phenos)

		plotdf.sort_values(by=['Phenotype', 'Gene set', 'Primary variant'], inplace=True)
		plotdf.reset_index(drop=True, inplace=True)

		# Add some variables for plotting
		gs_idx_map={}
		for g in gs:
			gs_idx_map[g]=gs.index(g)
		plotdf['variable_idx']=plotdf['Gene set'].map(gs_idx_map).astype(int)
		pheno_idx_map={}
		rev_pheno=quant_phenos.copy()
		rev_pheno.reverse()
		for p in quant_phenos:
			pheno_idx_map[p]=rev_pheno.index(p)
		plotdf['phenotype_idx']=plotdf.Phenotype.map(pheno_idx_map).astype(int)

		colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
		cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
		biggest_est=plotdf.Estimate.abs().max()
		norm=matplotlib.colors.Normalize(vmin=-biggest_est, vmax=biggest_est)
		plotdf['estimate_color']=plotdf.Estimate.apply(lambda x: cmap(norm(x)))

		plotdf['significance_color']='k'
		plotdf.loc[plotdf['BH FDR']<=0.05, 'significance_color']='red'

		plotdf['neglogp']=-np.log10(plotdf['p value'])
		plotdf['radius']=0.5*(plotdf.neglogp)/(plotdf.neglogp.max())

		plotdf['theta1']=plotdf['Primary variant'].map({'DBD Tier 1 SNVs':90, 'Large rare deletions':180, 'Large rare duplications':0, 'Any variant':270}).astype(int)
		plotdf['theta2']=plotdf.theta1+90

		# Make interaction plot
		fig, axs = plt.subplots(ncols=2, figsize=(10, 5), sharex=True, sharey=True)
		for idx, row in plotdf.iterrows():
			axs[0].add_artist(mpatches.Wedge((row['variable_idx'], row['phenotype_idx']), row['radius'], row['theta1'], row['theta2'], fc=row['estimate_color'], ec=row['significance_color'], lw=1))
		axs[0].set_xlim(-1, len(gs))
		axs[0].set_ylim(-1, 8)
		axs[0].set_aspect('equal', adjustable='box')

		axs[0].set_yticks([i for i in range(0, 8)], rev_pheno)
		axs[0].set_xticks([i for i in range(0, len(gs))], gs, rotation=90)

		# Make legend panel
		legdf=pd.DataFrame({'pvalue':[0.1, 0.01, 0.001, 0.0001, 0.0001, 0.001, 0.001, 0.001, 0.001], 'theta1':[90, 90, 90, 90, 90, 90, 180, 0, 270], 'y':[3, 2, 1, 2, 1, 7, 6.5, 5, 4.5], 'x':[0, 0, 0, 3, 3, 0, 0, 0, 0],
							'fillcolor':['#CCCCCC']*9, 'edgecolor':['k']*4+['red']+['k']*4,
							'text':['0.1', '0.01', '0.001', '0.0001', 'FDR<0.05', 'DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'Any variant']})
		legdf['neglogp']=-np.log10(legdf.pvalue)
		legdf['radius']=0.5*(legdf.neglogp)/(plotdf.neglogp.max())
		legdf['theta2']=legdf.theta1+90

		# Adjustments
		legdf.loc[legdf.theta1.isin([0, 270]), 'x']=legdf[legdf.theta1.isin([0, 270])].x-legdf[legdf.theta1.isin([0, 270])].radius

		legdf['text_x']=legdf.x+0.1
		legdf['text_y']=legdf.y+(0.5*legdf.radius)

		# Text adjustments
		legdf.loc[legdf.theta1.isin([180, 270]), 'text_y']=legdf[legdf.theta1.isin([180, 270])].y-(0.5*legdf[legdf.theta1.isin([180, 270])].radius)
		legdf.loc[legdf.theta1.isin([0, 270]), 'text_x']=legdf[legdf.theta1.isin([0, 270])].x+0.1+legdf[legdf.theta1.isin([0, 270])].radius

		for idx, row in legdf.iterrows():
			axs[1].add_artist(mpatches.Wedge((row['x'], row['y']), row['radius'], row['theta1'], row['theta2'], fc=row['fillcolor'], ec=row['edgecolor'], lw=1))
			# Add text
			axs[1].text(row['text_x'], row['text_y'], row['text'], color='k', va='center', ha='left')
		axs[1].set_aspect('equal', adjustable='box')

		# Colorbar
		scalarMap=cmx.ScalarMappable(norm=norm, cmap=cmap)
		scalarMap.set_array([])
		steps=[-biggest_est]
		for i in range(4):
			steps.append(steps[i]+(biggest_est/2.5))
		steps.append(biggest_est)
		plt.colorbar(scalarMap, ax=axs[1], ticks=steps, orientation='vertical')

		title='WGCNA modules'
		if gs==pathway:
			title='Pathways'
		plt.suptitle(title)
		plt.tight_layout()
		pdf.savefig()
		plt.close()
	pdf.close()
