import pandas as pd

# Parse annotated CNV calls

# Input and output
INPUT_DIR='../data/intracohort_gene_filter/Gene_anno'
OUTPUT_FILE='../data/filtered_calls/iPSC_cohort_calls.csv'

annos=['upstream', 'utr5', 'intron', 'exon', 'utr3', 'downstream']

# Load data
df=pd.DataFrame({'Sample':[], 'Chr':[], 'Start':[], 'End':[], 'CNV_type':[], 'Length':[], 'gnomAD_freq':[], 'Intracohort_count':[]})
for anno in annos:
	adf=pd.read_csv(f'{INPUT_DIR}/{anno}_annotated.bed', sep='\t', header=None,
					names=['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq', 'Intracohort_count',
							'Chr2', 'Start2', 'End2', 'Gene_ID', 'UNK', 'strand'])
	adf=adf[['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq', 'Intracohort_count', 'Gene_ID']]
	adf.columns=['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq', 'Intracohort_count', anno]
	adf=adf.groupby(['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq', 'Intracohort_count']).agg(lambda x: ';'.join(sorted(list(set(x)))))
	adf.reset_index(inplace=True)
	df=pd.merge(df, adf, on=['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Length', 'gnomAD_freq', 'Intracohort_count'], how='outer')

def cat_lst(x):
	exp_lst=[]
	for i in x.to_list():
		exp_lst+=i.split(';')
	return ';'.join(sorted(list(set([i for i in exp_lst if i not in ['.', '']]))))

# Merge by location
df.fillna('', inplace=True)
df=df.groupby(['Chr', 'Start', 'End', 'Sample', 'CNV_type', 'Intracohort_count']).agg({'Length':'first',
																						'gnomAD_freq':'max',
																						'upstream':cat_lst, 'utr5':cat_lst,
																						'intron':cat_lst, 'exon':cat_lst,
																						'utr3':cat_lst, 'downstream':cat_lst})
df.reset_index(inplace=True)

# Drop CNVs that do not intersect any part of genes
print(df.shape)
df=df[(~df[annos].isin(['.', ''])).any(axis=1)]
print(df.shape)

df.reset_index(inplace=True, drop=True)

# Annotate other samples with overlapping CNVs based on reciprocal overlap
# Check for 50% reciprocal overlap
samps=sorted(list(df.Sample.unique()))
outlst=[]
for idx, row in df.iterrows():
	if idx%100==0:
		print(idx)
	chr=row.Chr
	start=row.Start
	end=row.End
	CN=row.CNV_type
	length=row.Length

	# Find overlapping CNVs
	sdf=df[(df.Chr==chr) & (df.CNV_type==CN) & (((df.Start<=start) & (df.End>=start)) | ((df.Start<=end) & (df.End>=end)) | ((df.Start<=end) & (df.End>=start)))].copy()
	sdf['RS']=start
	sdf['RE']=end
	sdf['overlap']=(sdf[['End', 'RE']].min(axis=1))-(sdf[['Start', 'RS']].max(axis=1))

	sdf=sdf[sdf.overlap>=(0.5*sdf.Length)]
	sdf=sdf[sdf.overlap>=(0.5*length)]

	over_samps=sorted(list(sdf.Sample.unique()))

	out=[]
	for samp in samps:
		if samp in over_samps:
			out.append(sdf[sdf.Sample==samp].shape[0])
		else:
			out.append(0)
	outlst.append(out)
outdf=pd.DataFrame(outlst, columns=samps)
df=pd.merge(df, outdf, right_index=True, left_index=True)

# Filter for intracohort frequency <= 5
df['Intracohort_man']=(df[samps]>0).sum(axis=1)
df=df[df.Intracohort_man<=5]
print(df.shape)

# Save
df.to_csv(OUTPUT_FILE, index=False)

print(df.gnomAD_freq.max())
print(df[df.gnomAD_freq<=0.001].shape)
