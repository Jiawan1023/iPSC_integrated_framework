import pandas as pd
import numpy as np

# Parse STR calls and filter for:
# 1. Calls that are rare in Ensembl (<1% frequency)
# 2. Calls that are not rare, but represent expansions >2SD than the Ensembl cohort mean

# Input and output
INPUT_DIR="../data/GangSTR/"
ENSEMBL_DIR="/data7/johnathan/ensembleTR/statstr"
OUTPUT_DIR="../data/Filtered_STRs/"

# Load data
df=pd.DataFrame()
for chrom in range(1, 23):
	cdf=pd.read_csv(f'{INPUT_DIR}/chr{chrom}/chr{chrom}_GT.tsv', sep='\t')
	df=pd.concat([df, cdf])

# Parse period length
per_col=[i for i in df.columns.to_list() if 'PERIOD' in i][0]
df['Period']=df[per_col].str.split(',', expand=True)[0].astype(int)
del df[per_col]

# Rename columns
new_cols=[i.split(']')[1].split(':')[0] if i!='Period' else i for i in df.columns.to_list()]
df.columns=new_cols

df['Locus']=df.CHROM+':'+df.POS.astype(str)

# Load ENSEMBL data
edf=pd.DataFrame()
for chrom in range(1, 23):
	cdf=pd.read_csv(f'{ENSEMBL_DIR}/ensemble_chr{chrom}_filtered.tab', sep='\t')
	# Restrict to sites in data
	cdf['Locus']=cdf.chrom+':'+cdf.start.astype(str)
	cdf=cdf[cdf.Locus.isin(df.Locus.to_list())]
	edf=pd.concat([edf, cdf])

# Some loci are duplictaed - drop duplicates
edf=edf[~edf.duplicated(subset=['Locus'], keep=False)]

df=df[df.Locus.isin(edf.Locus.to_list())]

# Split ENSEMBL data by allele
ensembl_stats=edf[['Locus', 'mean', 'var']]
ensembl_stats['SD']=ensembl_stats['var']**0.5
ensembl_stats=ensembl_stats[['Locus', 'mean', 'SD']]
ensembl_stats.columns=['Locus', 'ENSEMBL_Mean', 'ENSEMBL_SD']

edf=edf[['Locus', 'afreq']]
edf.afreq=edf.afreq.str.split(',')
edf=edf.explode('afreq')
edf[['ALT', 'ENSEMBL_freq']]=edf.afreq.str.split(':', expand=True)
edf.ALT=edf.ALT.str.upper()
edf.ENSEMBL_freq=edf.ENSEMBL_freq.astype(float)
edf=edf[['Locus', 'ALT', 'ENSEMBL_freq']]

# Explode by sample
df=df.melt(id_vars=['CHROM', 'POS', 'REF', 'ALT', 'Period', 'Locus'], var_name='Sample', value_name='GT')

# Drop REF alleles
df=df[df.GT!='0/0']
df=df[df.GT!='.']

# Annotate genotype info
df['N_alleles']=df.GT.apply(lambda x: len(list(set([i for i in x.split('/') if i!='0']))))
df['Zygosity']='COMPOUND_HET'
df.loc[df.GT.str.contains('0'), 'Zygosity']='HET'
df.loc[(df.N_alleles==1) & (~df.GT.str.contains('0')), 'Zygosity']='HOM'

# Explode by allele
df.GT=df.GT.str.split('/')
df=df.explode('GT')
df=df[df.GT!='0']
df.GT=df.GT.astype(int)
df.drop_duplicates(inplace=True)

df=df[['CHROM', 'POS', 'Locus', 'REF', 'ALT', 'Period', 'Sample', 'GT', 'Zygosity']]

# Extract correct ALT allele
df['ALT']=df[['ALT', 'GT']].apply(lambda row: row['ALT'].split(',')[row.GT-1], axis=1)

# Get length of ref and alt alleles
df['REF_LEN']=(df.REF.str.len()/df.Period).astype(int)
df['ALT_LEN']=(df.ALT.str.len()/df.Period).astype(int)

# Convert REF and ALT to upper
df.REF=df.REF.str.upper()
df.ALT=df.ALT.str.upper()

# Merge with ENSEMBL data
df=pd.merge(df, edf, on=['Locus', 'ALT'], how='left')
df.fillna(-1, inplace=True)
df=pd.merge(df, ensembl_stats, on=['Locus'], how='left')

# Annotate frequency in the cohort
str_counts=pd.DataFrame(df[['Locus', 'ALT', 'Sample']].groupby(['Locus', 'ALT']).agg(lambda x: len(list(set(x)))))
str_counts.columns=['Intracohort_count']
str_counts.reset_index(inplace=True)

df=pd.merge(df, str_counts, on=['Locus', 'ALT'], how='left')

# Annotate with filters
df['Rare_ENSEMBL']=((df.ENSEMBL_freq<=0.01) & (df.ENSEMBL_freq>-1))
df['Expansion_ENSEMBL']=(df.ALT_LEN>=(df.ENSEMBL_Mean+(2*df.ENSEMBL_SD)))
df['Cohort_freq']=(df.Intracohort_count<=5)

# Save unfiltered data
df.to_csv(f'{OUTPUT_DIR}/Unfiltered_loci.csv', index=False)

# Filter data
print(df.shape)
df=df[(df.Cohort_freq)]
print(df.shape)
df=df[(df.Rare_ENSEMBL) | (df.Expansion_ENSEMBL)]
print(df.shape)

df=df[['Sample', 'Locus', 'REF', 'ALT', 'REF_LEN', 'ALT_LEN', 'Zygosity', 'ENSEMBL_freq', 'ENSEMBL_Mean', 'ENSEMBL_SD', 'Intracohort_count']]

df.REF=df.REF.str.lower()
df.ALT=df.ALT.str.lower()

df.to_csv(f'{OUTPUT_DIR}/Filtered_loci.csv', index=False)
