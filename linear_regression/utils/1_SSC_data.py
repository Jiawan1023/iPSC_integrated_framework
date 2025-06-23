import pandas as pd

import os
cwd = os.getcwd()
dropbox=cwd.split('Dropbox')[0]+'Dropbox/'
dropbox=dropbox.replace('\\', '/')

# Load in SSC head circumference data
ssc_hwhc=pd.read_csv(dropbox+"SSC data/Phenotype Data Set 9/Proband Data/ssc_hwhc.csv")
ssc_cd=pd.read_csv(dropbox+"SSC data/Phenotype Data Set 9/Proband Data/ssc_core_descriptive.csv")

ssc=pd.merge(ssc_hwhc[['individual', 'head_circumference']], ssc_cd[['individual', 'sex', 'age_at_ados']], how='inner')
ssc=ssc[~ssc.isnull().any(axis=1)]

# Convert head circumference data to Z-score
mhc=pd.read_csv('files/Population_Data/Male_HC.csv', index_col=0)
fhc=pd.read_csv('files/Population_Data/Female_HC.csv', index_col=0)

mhc['Mean']=mhc['50th']
fhc['Mean']=fhc['50th']

mhc['SD']=(((mhc['3rd']-mhc['Mean'])/-1.881)+((mhc['10th']-mhc['Mean'])/-1.282)+((mhc['25th']-mhc['Mean'])/-0.674)+((mhc['75th']-mhc['Mean'])/0.674)+((mhc['90th']-mhc['Mean'])/1.282)+((mhc['97th']-mhc['Mean'])/1.881))/6
fhc['SD']=(((fhc['3rd']-fhc['Mean'])/-1.881)+((fhc['10th']-fhc['Mean'])/-1.282)+((fhc['25th']-fhc['Mean'])/-0.674)+((fhc['75th']-fhc['Mean'])/0.674)+((fhc['90th']-fhc['Mean'])/1.282)+((fhc['97th']-fhc['Mean'])/1.881))/6

# Round Age to nearest in HC data
ssc['nearest_age']=ssc.age_at_ados.apply(lambda x: min(mhc.index.to_list(), key=lambda y:abs(y-x)))
ssc['Mean']=ssc.nearest_age.map(mhc.Mean.to_dict())
ssc['SD']=ssc.nearest_age.map(mhc.SD.to_dict())
ssc.loc[ssc.sex=='female', 'Mean']=ssc[ssc.sex=='female'].nearest_age.map(fhc.Mean.to_dict())
ssc.loc[ssc.sex=='female', 'SD']=ssc[ssc.sex=='female'].nearest_age.map(fhc.SD.to_dict())

ssc['HC_zscore']=(ssc.head_circumference-ssc.Mean)/ssc.SD

# Load in other data from WGS paper
ssc_wgs=pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/CLEAN_CODE/SSC/SSC_Data/SSC_whole_cohort.csv")

ssc=pd.merge(ssc_wgs, ssc, right_on='individual', left_on='Sample', how='left')

# Note which kinds of variants are available
ssc['SNV']=True
ssc.loc[ssc['Missense'].isnull(), 'SNV']=False
ssc['CNV']=True
ssc.loc[ssc['Genes del.'].isnull(), 'CNV']=False
ssc['STR']=True
ssc.loc[ssc.STRs.isnull(), 'STR']=False

ssc=ssc[['Sample', 'DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'Any variant', 'SNV', 'CNV', 'STR', 'Sex', 'Age',
			'Full scale IQ', 'Externalizing behavior (ABCL/CBCL)', 'Internalizing behavior (ABCL/CBCL)', 'Social responsiveness (SRS)', 'Repetitive behavior (RBS-R)', 'Coordination disorder (DCDQ)', 'BMI z-score', 'HC_zscore',
			'head_circumference', 'nearest_age', 'Mean', 'SD']]
ssc.sort_values(by='Sample', inplace=True)

# Save
ssc.to_csv('files/SSC_phenotypes.csv', index=False)
