import pandas as pd
import numpy as np

import scipy.stats as stats
import statsmodels.formula.api as smf

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams['pdf.fonttype']=42

import re

# Use linear regression to test for the effects of module/pathway genes, a first hit, and the interaction of module/pathway genes and a first hit

# Load in data
svip_pheno=pd.read_csv('files/Searchlight_phenotypes.csv')
svip_gs=pd.read_csv('files/Searchlight_gene_set_counts.csv')

svip=pd.merge(svip_pheno, svip_gs, on='Sample', how='inner')

# Based on the distribution of counts, treat WGCNA module gene sets as quantitative and pathway gene sets as binary
wgcna=[i for i in svip_gs.columns.to_list() if '_' not in i and i!='Sample']
pathway=[i for i in svip_gs.columns.to_list() if '_' in i]

phenos=['Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Autism behavior (BSI)', 'BMI z-score', 'Head circumference z-score']

replace_string='[ .()/-]'
def run_model(moddf, output_col, cohort, gene_set):
	# Replace problem characters in variable names
	form='%s ~ Age + Sex + %s' % (re.sub(replace_string, '_', output_col), gene_set)
	mod=smf.ols(formula=form, data=moddf)
    
	res=mod.fit()
	
	# Parse model
	ci=res.conf_int(alpha=0.05)
	num_vars=ci.shape[0]
	ci.columns=['95% C.I. lower', '95% C.I. upper']
	
	r2=res.rsquared
	
	mod_name='%s ~ Age + Sex + %s' % (output_col, gene_set)
	
	res_dict={'Primary variant':[cohort]*num_vars, 'Gene set':[gene_set]*num_vars, 'Phenotype':[output_col]*num_vars, 'Model':[mod_name]*num_vars,
				'Variable':['Intercept', 'Age', 'Sex', gene_set], 'R2':[r2]*num_vars}
	mod_res=pd.DataFrame(res_dict)
	mod_res.index=['Intercept', 'Age', 'Sex', gene_set]
	
	mod_res['Estimate']=res.params
	mod_res['Error']=res.bse
	mod_res['p value']=res.pvalues
	mod_res=pd.merge(mod_res, ci, left_index=True, right_index=True)
	
	mod_res=mod_res[['Primary variant', 'Gene set', 'Phenotype', 'Model', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]
	return mod_res

outdf=pd.DataFrame()
for pheno in phenos:
	for cohort in ['16p11.2 deletion', '16p11.2 duplication']:
		for gs in wgcna+pathway:
			moddf=svip[svip.Cohort==cohort][[pheno, gs, 'Age', 'Sex']].copy()
			moddf=moddf[~moddf.isnull().any(axis=1)]
			
			# If the number of samples with a variant in a gene set is very small - skip
			if moddf[moddf[gs]>0].shape[0]<2:
				print('Skipping', pheno, cohort, gs, '- too few samples')
				continue

			# Clean up binary variables
			moddf['Sex']=moddf.Sex.map({'M':0, 'F':1})
			if gs in pathway:
				moddf.loc[moddf[gs]>0, gs]=1
			
			# Scale quantiative variables
			for qv in [pheno, 'Age', gs]:
				if qv==gs and gs not in wgcna:
					continue
				moddf[qv]=(moddf[qv]-moddf[qv].mean())/moddf[qv].std()
			
			# Replace problem characters in variable names
			cols=moddf.columns.to_list()
			moddf.columns=[re.sub(replace_string, '_', i) for i in cols]
			
			out=run_model(moddf, pheno, cohort, gs)
			out['Sample size']=moddf.shape[0]
			out.reset_index(inplace=True, drop=True)
			outdf=pd.concat([outdf, out])
outdf=outdf[['Primary variant', 'Gene set', 'Phenotype', 'Model', 'Sample size', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]

# FDR correction
outdf['BH FDR']=stats.false_discovery_control(outdf['p value'], method='bh')
print(outdf)

# Save
outdf.to_csv('Result_tables/12_Searchlight_regression.csv', index=False)

# Make figures
pdf=PdfPages('Figures/12_Searchlight_results.pdf')
i=0
for cohort in ['16p11.2 deletion', '16p11.2 duplication']:
	sup_idx=i
	sup_title=cohort
	fig, axs=plt.subplots(ncols=2, figsize=(10, 5))
	for gs in [wgcna, pathway]:
		idx=[wgcna, pathway].index(gs)
		
		plotdf=outdf[(outdf['Primary variant']==cohort) & (outdf.Variable.isin(gs))].copy()
		plotdf['star']=''
		plotdf.loc[plotdf['p value']<=0.05, 'star']='*'
		plotdf.loc[plotdf['BH FDR']<=0.05, 'star']='**'
		
		heatdf=plotdf.pivot(index='Gene set', columns='Phenotype', values='Estimate')
		stardf=plotdf.pivot(index='Gene set', columns='Phenotype', values='star')
		
		gidx=[i for i in gs if i in plotdf['Gene set'].to_list()]
		
		heatdf=heatdf.loc[gidx, phenos]
		stardf=stardf.loc[gidx, phenos]
		
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
	
	i+=1
pdf.close()

