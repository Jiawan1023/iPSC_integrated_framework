import pandas as pd
import numpy as np

import scipy.stats as stats
import statsmodels.formula.api as smf

import re

# Use linear regression to test for the effects of module/pathway genes, a first hit, and the interaction of module/pathway genes and a first hit

# Load in data
pheno=pd.read_csv('files/16p12_phenotypes.csv')
gs=pd.read_csv('files/16p12_gene_set_counts.csv')

df=pd.merge(pheno, gs, on='Sample', how='inner')

# Rename age column
df['Age']=df['Age (years)']

# Based on the distribution of counts, treat WGCNA module gene sets as quantitative and pathway gene sets as binary
wgcna=['module_gene_name']
pathway=['pathway_gene_name']

quant_phenos=['De Vries Score', 'BMI Z Score', 'Head Circumference Z Score', 'HRS-MAT', 'SRS Raw Score']

replace_string='[ .()/-]'
def run_model(moddf, output_col, gene_set, regtype='Linear'):
	# Replace problem characters in variable names
	form='%s ~ Age + Sex + %s' % (re.sub(replace_string, '_', output_col), gene_set)
	if regtype=='Linear':
		mod=smf.ols(formula=form, data=moddf)
	else:
		mod=smf.logit(formula=form, data=moddf)
    
	res=mod.fit()
	
	# Parse model
	ci=res.conf_int(alpha=0.05)
	num_vars=ci.shape[0]
	ci.columns=['95% C.I. lower', '95% C.I. upper']
	
	if regtype=='Linear':
		r2=res.rsquared
	else:
		r2=res.prsquared
	
	mod_name='%s ~ Age + Sex + %s' % (output_col, gene_set)
	
	res_dict={'Gene set':[gene_set]*num_vars, 'Phenotype':[output_col]*num_vars, 'Regression Type':[regtype]*num_vars, 'Model':[mod_name]*num_vars,
				'Variable':['Intercept', 'Age', 'Sex', gene_set], 'R2':[r2]*num_vars}
	mod_res=pd.DataFrame(res_dict)
	mod_res.index=['Intercept', 'Age', 'Sex', gene_set]
	
	mod_res['Estimate']=res.params
	mod_res['Error']=res.bse
	mod_res['p value']=res.pvalues
	mod_res=pd.merge(mod_res, ci, left_index=True, right_index=True)
	
	mod_res=mod_res[['Gene set', 'Phenotype', 'Regression Type', 'Model', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]
	return mod_res

outdf=pd.DataFrame()
for pheno in quant_phenos:
	for gs in wgcna+pathway:
		moddf=df[[pheno, gs, 'Age', 'Sex']].copy()
		moddf=moddf[~moddf.isnull().any(axis=1)]
		
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
		
		regtype='Linear'
		
		out=run_model(moddf, pheno, gs, regtype)
		out['Sample size']=moddf.shape[0]
		out.reset_index(inplace=True, drop=True)
		outdf=pd.concat([outdf, out])
outdf=outdf[['Gene set', 'Phenotype', 'Regression Type', 'Model', 'Sample size', 'Variable', 'Estimate', 'Error', '95% C.I. lower', '95% C.I. upper', 'p value', 'R2']]
print(outdf)

# FDR correction
outdf['BH FDR']=stats.false_discovery_control(outdf['p value'], method='bh')
print(outdf)

print(outdf[(outdf['p value']<=0.05) & (outdf.Variable.isin(wgcna+pathway))])

# Save
outdf.to_csv('Result_tables/9_16p12_regression.csv', index=False)
