import pandas as pd
import numpy as np

import scipy.stats as stats
import statsmodels.formula.api as smf

import re

# Use linear regression to test for the effects of module/pathway genes, a first hit, and the interaction of module/pathway genes and a first hit

# Load in data
ssc_pheno=pd.read_csv('files/SSC_phenotypes.csv')
ssc_gs=pd.read_csv('files/SSC_gene_set_counts.csv')

ssc=pd.merge(ssc_pheno, ssc_gs, on='Sample', how='inner')

# Based on the distribution of counts, treat WGCNA module gene sets as quantitative and pathway gene sets as binary
wgcna=['#input_module_names']
pathway=['input_pathway_names']

phenos=['Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score', 'HC_zscore']

replace_string='[ .()/-]'
def run_model(moddf, output_col, primary_variant, gene_set):
	# Replace problem characters in variable names
	form='%s ~ Age + Sex + %s + %s + %s*%s' % (re.sub(replace_string, '_', output_col), re.sub(replace_string, '_', primary_variant), gene_set,
												re.sub(replace_string, '_', primary_variant), gene_set)
	mod=smf.ols(formula=form, data=moddf)
    
	res=mod.fit()
	
	# Parse model
	ci=res.conf_int(alpha=0.05)
	num_vars=ci.shape[0]
	ci.columns=['95% C.I. lower', '95% C.I. upper']
	
	r2=res.rsquared
	
	mod_name='%s ~ Age + Sex + %s + %s + %s*%s' % (output_col, primary_variant, gene_set, primary_variant, gene_set)
	
	res_dict={'Primary variant':[primary_variant]*num_vars, 'Gene set':[gene_set]*num_vars, 'Phenotype':[output_col]*num_vars, 'Model':[mod_name]*num_vars,
				'Variable':['Intercept', 'Age', 'Sex', primary_variant, gene_set, 'Interaction'], 'R2':[r2]*num_vars}
	mod_res=pd.DataFrame(res_dict)
	mod_res.index=['Intercept', 'Age', 'Sex', re.sub(replace_string, '_', primary_variant), gene_set, re.sub(replace_string, '_', primary_variant)+':'+gene_set]
	
	mod_res['Estimate']=res.params
	mod_res['Error']=res.bse
	mod_res['p value']=res.pvalues
	mod_res=pd.merge(mod_res, ci, left_index=True, right_index=True)
	
	mod_res=mod_res[['Primary variant', 'Gene set', 'Phenotype', 'Model', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]
	return mod_res

outdf=pd.DataFrame()
for pheno in phenos:
	for fh in ['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'Any variant']:
		for gs in wgcna+pathway:
			moddf=ssc[[pheno, gs, fh, 'Age', 'Sex']].copy()
			moddf=moddf[~moddf.isnull().any(axis=1)]
			
			# Clean up binary variables
			moddf[fh]=moddf[fh].astype(int)
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
			
			out=run_model(moddf, pheno, fh, gs)
			out['Sample size']=moddf.shape[0]
			out.reset_index(inplace=True, drop=True)
			outdf=pd.concat([outdf, out])
outdf=outdf[['Primary variant', 'Gene set', 'Phenotype', 'Model', 'Sample size', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]
print(outdf)

# FDR correction
outdf['BH FDR']=stats.false_discovery_control(outdf['p value'], method='bh')
print(outdf)


# Save
outdf.to_csv('Result_tables/4_SSC_regression.csv', index=False)
