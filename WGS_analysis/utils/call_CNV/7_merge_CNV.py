import pandas as pd

# Merge and QC filter CNVs

# Filter
# 1. Remove sex chromosome CNVs
# 2. Remove CNVs with q0>0.5 - See: https://github.com/abyzovlab/CNVpytor/issues/115
# 3. Remove CNVs with pN>0.5
# 3. Remove CNVs less than 3 bins long (1500bp)
# 4. Remove CNVs with e_val1 > 0.00001

# Merge CNVs that are:
# 1. Separated by less than 5kb, and
# 2. Less than 20% of combined CNV length

# Input and output
INPUT_CNV='../data/region_filter/UCSC_region_filter/GRC_Exclusions_filter.bed'
OUTPUT_DIR='../data/CNV_merge'

# Load data
df=pd.read_csv(INPUT_CNV, sep='\t', header=None, names=['Chr', 'Start', 'End', 'Length', 'CNV_type', 'read_depth', 'e_val1',
														'e_val2', 'e_val3', 'e_val4', 'q0', 'pN', 'dG', 'gene_ids', 'gene_names',
														'Sample', 'gnomAD_freq'])

# Filter
print(df.shape)
df=df[df.Chr!='chrY']
print(df.shape)
df=df[df.q0<0.5]
print(df.shape)
df=df[df.pN<0.5]
print(df.shape)
df=df[df.Length>=1500]
print(df.shape)
df=df[df.e_val1<=0.00001]
print(df.shape)

# Collapse duplicated annotations
df=df.groupby(['Chr', 'Start', 'End', 'Sample', 'CNV_type']).agg({'Length':'first',
																	'read_depth':'first',
																	'e_val1':'first',
																	'e_val2':'first',
																	'e_val3':'first',
																	'e_val4':'first',
																	'q0':'max',
																	'pN':'max',
																	'dG':'first',
																	'gnomAD_freq':'max'})
df.reset_index(inplace=True)
print(df.shape)

# Save filtered calls pre-merged
df.to_csv(f'{OUTPUT_DIR}/QC_filtered_CNVs.csv', index=False)

# Merge CNVs
df.sort_values(by=['Sample', 'CNV_type', 'Chr', 'Start', 'End'], inplace=True)
df.reset_index(inplace=True, drop=True)

gap=5000
while True:
	# Shift data
	df['Chr_shift']=df.Chr.shift(-1)
	df['Sample_shift']=df.Sample.shift(-1)
	df['CN_shift']=df.CNV_type.shift(-1)
	df['start_shift']=df.Start.shift(-1)
	df['end_shift']=df.End.shift(-1)

	df['start_min']=df[['Start', 'start_shift']].min(axis=1)
	df['end_max']=df[['End', 'end_shift']].max(axis=1)

	# Define overlaps
	df['overlap']=(((df.start_shift-gap)<=df.End) & (df.Sample==df.Sample_shift) & (df.Chr==df.Chr_shift) & (df.CNV_type==df.CN_shift) & ((df.start_shift-df.End)<=0.2*(df.end_max-df.start_min)))
	df['is_overlapped']=df.overlap.shift(1)
	df.loc[df.is_overlapped.isnull(), 'is_overlapped']=False

	# Check if overlaps still exist
	# If not, exit
	if df[df.overlap].shape[0]==0:
		break

	# Resolve overlaps
	df.loc[df.overlap, 'Start']=df[df.overlap].start_min
	df.loc[df.overlap, 'End']=df[df.overlap].end_max

	# Drop unneeded columns
	df=df[~((df.is_overlapped) & (~df.overlap))]

df['Length']=df.End+1-df.Start

# Drop unneeded columns
df=df[['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq']]

# Save as BED
for cn in ['del', 'dup']:
	df[df.CNV_type.str.contains(cn)].to_csv(f'{OUTPUT_DIR}/{cn}_merged.bed', sep='\t', header=False, index=False)
