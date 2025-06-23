import pandas as pd

import os
cwd = os.getcwd()
dropbox=cwd.split('Dropbox')[0]+'Dropbox/'
dropbox=dropbox.replace('\\', '/')

# Load in 16p12 proband data
wgs=pd.read_csv(dropbox+"16p12.2 project/Human patients project/WGS paper/CLEAN_CODE/DD_cohort/Table_S1A.csv")

wgs=wgs[(wgs['Estonian Biobank Sample']!='X') & (wgs['16p12.1 deletion']=='Carrier') & (wgs.WGS=='X')]
wgs=wgs[['Sample', 'Family', 'Relationship', 'Age (years)', 'Sex', 'WGS', 'Microarray', 'Phenotypic Domains',
		'Behavioral features (Child domain)', 'Psychiatric features (Child domain)', 'Nervous System Abnormalities (Child domain)', 'Congenital Anomalies (Child domain)', 'Growth/Skeletal Defects (Child domain)',
		'De Vries Score', 'BMI Z Score', 'Head Circumference Z Score', 'HRS-MAT', 'SRS Raw Score']]

# Drop samples without any phenotype information
wgs=wgs[~wgs[bin_phenos+quant_phenos].isnull().all(axis=1)]

# Remove samples without age or sex information
wgs=wgs[~wgs[['Age (years)', 'Sex']].isnull().any(axis=1)]
# Remove stillborn sample
wgs=wgs[wgs['Age (years)']!='Stillborn']

# Update sample IDs to SG codes
df=pd.read_excel(dropbox+"16p12.2 project/Human patients project/WGS paper/11_Variant Integration/Supp Table 1/Extra_Sample_Information.xlsx")
df.index=df.New_Sample_ID
wgs['Sample']=wgs.Sample.map(df.Sample.to_dict())

# Save
wgs.to_csv('files/16p12_phenotypes.csv', index=False)
