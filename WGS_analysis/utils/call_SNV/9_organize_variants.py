import pandas as pd
import networkx as nx

# Clean up variants and organize

chrs=[str(i) for i in range(1, 23)]+['X', 'Y']

df=pd.DataFrame()
for i in chrs:
        chrdf=pd.read_csv('../data/HAIL/tables/chr'+i+'.csv')
        df=pd.concat([df, chrdf])
print(df.shape)

# Remove any variants present in unrelated individuals
def pedgraph(pedfile):
	ped=pd.read_csv(pedfile, sep='\t', header=None, names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])

	# Create pedigrees as a graph
	pedg=nx.DiGraph()
	pedg.add_nodes_from(ped.IID.to_list())

	# Add parent-child edges
	parent=[tuple(r) for r in ped[ped.Father!='0'][['Father', 'IID']].to_numpy()]+[tuple(r) for r in ped[ped.Mother!='0'][['Mother', 'IID']].to_numpy()]

	pedg.add_edges_from(parent)
	return pedg

def find_rels(samp, nxgraph, to_rem):
	preds=[samp]
	# Find all predecessors of current sample
	while True:
		old_preds=preds.copy()
		for p in old_preds:
			preds+=list(nxgraph.predecessors(p))
		preds=sorted(list(set(preds)))
		if len(preds)==len(old_preds):
			break

	# Find all successors of predecessors
	succ=[]
	for p in preds:
		succ+=list(nxgraph.successors(p))
	while True:
		old_succ=succ.copy()
		for s in old_succ:
			succ+=list(nxgraph.successors(s))
		succ=sorted(list(set(succ)))
		if len(succ)==len(old_succ):
			break

	rels=sorted(list(set(preds+succ)))

	# Some sample IDs need to be removed because they are placeholders to make the correct graph shape
	rels=[i for i in rels if i not in to_rem]
	return rels

pedg=pedgraph('files/iPSC.ped') # Note that PED has several placeholder individuals to maintain relationships
sample_paths={}
for n in pedg.nodes:
	sample_paths[n]=find_rels(n, pedg, ['CR1F', 'CR1M', 'SMF', 'SMM'])

df['vid']=df.Chr+'_'+df.Pos.astype(str)+'_'+df.Ref+'_'+df.Alt
fdf=df[['Sample', 'vid']].drop_duplicates().groupby('vid')['Sample'].apply(list).reset_index(name='samples')
fdf.samples=fdf.samples.apply(lambda x: sorted(x))

def family_unique(samps):
    for s in samps:
        srels=sample_paths[s]
        if set(srels)>=set(samps):
            continue
        else:
            return False
    return True

fdf['family_unique']=fdf.samples.apply(family_unique)
fdf=fdf[fdf.family_unique]

df=df[df.vid.isin(fdf.vid.to_list())]
print(df.shape)

# Organize variants in 2 ways:
# 1. A burden table counting the number of variants of each type
# 2. A list of the variant-gene/transcripts (separated by coding/non-coding)

# Burden Table
def get_burden(df):
	exdf=df[['Sample', 'vid', 'variant_type', 'Chr']].copy()
	exdf.variant_type=exdf.variant_type.str.split(';')
	exdf=exdf.explode('variant_type')
	burden_df=exdf.groupby(['Sample', 'vid', 'variant_type']).apply(lambda x: 1, include_groups=False)
	burden_df=burden_df.reset_index()
	burden_df=burden_df[['Sample', 'variant_type', 0]].groupby(['Sample', 'variant_type']).sum()
	burden_df=burden_df.reset_index()
	burden_df.columns=['Sample', 'variant_type', 'burden']
	# Convert to wide
	burden_wide=pd.pivot(burden_df, index='Sample', columns='variant_type', values='burden')
	burden_wide.fillna(0, inplace=True)
	burden_wide=burden_wide.astype(int)
	# Because the total variant count may not equal the individual variant counts (single variant may affect multiple transcripts/have multiple effects), add in the total burden
	tot_burd=df[['Sample', 'vid', 'Chr']].groupby(['Sample', 'vid']).apply(lambda x: 1, include_groups=False)
	tot_burd=tot_burd.reset_index()
	tot_burd.columns=['Sample', 'vid', 'total_burden']
	tot_burd=tot_burd[['Sample', 'total_burden']].groupby(['Sample']).sum()
	burden_wide=pd.merge(burden_wide, tot_burd, left_index=True, right_index=True)
	# Clean up and save
	burden_wide=burden_wide[['total_burden', 'missense', 'lof', 'splice_lof', 'splice', 'utr3', 'utr5', 'upstream', 'downstream', 'intron']]
	burden_wide.reset_index(inplace=True)
	return burden_wide
# Save burden
burden_wide=get_burden(df)
burden_wide.to_csv('../data/variant_tables/all_burden_table.csv', index=False)
# Repeat, but only for coding transcripts
burden_wide=get_burden(df[df.biotype=='protein_coding'])
burden_wide.to_csv('../data/variant_tables/coding_burden_table.csv', index=False)

# Separate variants into separate tables for protein-coding and non-coding transcripts and coding/noncoding mutations
# For now, only use upstream/UTR noncoding variants
print(df.variant_type.value_counts())
variant_types=sorted(list(df.variant_type.unique()))
coding_muts=[i for i in variant_types if 'lof' in i or 'splice' in i or 'splice_lof' in i or 'missense' in i]
noncoding_muts=[i for i in variant_types if 'upstream' in i or 'downstream' in i or 'utr3' in i or 'utr5' in i or 'intron' in i]

df[(df.variant_type.isin(coding_muts)) & (df.biotype=='protein_coding')].to_csv('../data/variant_tables/transcript_tables/coding_mut_coding_transcript.csv', index=False)
df[(df.variant_type.isin(noncoding_muts)) & (df.biotype=='protein_coding')].to_csv('../data/variant_tables/transcript_tables/noncoding_mut_coding_transcript.csv', index=False)
df[(df.variant_type.isin(coding_muts)) & (df.biotype!='protein_coding')].to_csv('../data/variant_tables/transcript_tables/coding_mut_noncoding_transcript.csv', index=False)
df[(df.variant_type.isin(noncoding_muts)) & (df.biotype!='protein_coding')].to_csv('../data/variant_tables/transcript_tables/noncoding_mut_noncoding_transcript.csv', index=False)

# Condense by gene
df=df[['Sample', 'Chr', 'Pos', 'Ref', 'Alt', 'gene', 'gene_id', 'transcript', 'biotype',
		'variant_type', 'consequence', 'gnomad_freq', 'cohort_freq',
		'vid']].groupby(['Sample', 'Chr', 'Pos', 'Ref', 'Alt', 'gene', 'gene_id', 'gnomad_freq', 'cohort_freq', 'vid']).agg(lambda x: '|'.join(sorted(list(set([str(i) for i in x])))))
df.reset_index(inplace=True)
print(df)

# Split by coding and non-coding transcripts
df[df.biotype.str.contains('protein_coding[^_]', regex=True)].to_csv('../data/variant_tables/gene_tables/coding_gene_variants.csv', index=False)
df[~df.biotype.str.contains('protein_coding[^_]', regex=True)].to_csv('../data/variant_tables/gene_tables/noncoding_gene_variants.csv', index=False)
