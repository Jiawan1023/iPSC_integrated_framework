import pandas as pd

# Annotate STR variants with ANNOVAR annotations and filter for desired STR types

# Input and output
FILTERED_STR='../data/Filtered_STRs/Filtered_loci.csv'
ANNOVAR_ANNO='../data/ANNOVAR/Parsed_ANNOVAR.txt'
GTF='/data5/bx_reference/hg38/annotations/gene_annotations/GENCODE39/gencode.v39.annotation.gtf.gz'
OUTPUT_DIR='../data/Final_STR'

FUNC_LIST=['exonic', 'splicing', 'intronic', 'UTR3', 'UTR5', 'upstream', 'downstream']

# Load data
df=pd.read_csv(FILTERED_STR)
annodf=pd.read_csv(ANNOVAR_ANNO, sep='\t')

gtf=pd.read_csv(GTF, compression='gzip', sep='\t', comment='#', header=None,
					names=['Chr', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attributes'])
gtf=gtf[gtf.Feature=='gene']

# Rename annotation columns
annodf.columns=['CHROM', 'POS', 'REF', 'ALT', 'Func', 'Gene', 'ExonicFunc']

# Fix ANNOVAR characters
for col in ['Func', 'Gene', 'ExonicFunc']:
	annodf[col]=annodf[col].str.replace('\\x3b', ';')

# Remove unneeded annotations
annodf['keep']=False
for fun in FUNC_LIST:
	annodf.loc[annodf.Func.str.contains(fun), 'keep']=True
annodf=annodf[annodf.keep]

annodf=annodf[annodf.Func!='ncRNA_exonic']
annodf=annodf[annodf.Func!='ncRNA_splicing']
annodf=annodf[annodf.Func!='ncRNA_intronic']

# Merge annotations
annodf['Locus']=annodf.CHROM+':'+annodf.POS.astype(str)
print(df.shape)
df=pd.merge(df, annodf[['Locus', 'REF', 'ALT', 'Func', 'Gene', 'ExonicFunc']], on=['Locus', 'REF', 'ALT'], how='left')

# Drop any STRs that do not interact with genes
df=df[~df.Func.isnull()]
print(df.shape)

# Split single function, multigene sites
df.Gene=df[['Func', 'Gene']].apply(lambda row: [row.Gene] if ';' in row.Func else row.Gene.split(';'), axis=1)
df=df.explode('Gene')

# Make VID
df['Period']=(df.REF.str.len()/df.REF_LEN).astype(int)
df['RU']=df[['REF', 'Period']].apply(lambda row: row.REF[0:row.Period], axis=1).str.upper()

df['VID']=df.Locus+'_'+df.REF_LEN.astype(str)+df.RU+'_'+df.ALT_LEN.astype(str)+df.RU
df.VID=df.VID.str.replace(':', '_')

df=df[['Sample', 'VID', 'REF_LEN', 'ALT_LEN', 'Zygosity', 'ENSEMBL_freq', 'ENSEMBL_Mean', 'ENSEMBL_SD', 'Intracohort_count',
	   	'Func', 'Gene', 'ExonicFunc']]

# Annotate gene types
gtf['Gene_name']=gtf.Attributes.str.split('gene_name ', expand=True)[1].str.split(';', expand=True)[0].str.replace('"', '')
gtf['Gene_type']=gtf.Attributes.str.split('gene_type ', expand=True)[1].str.split(';', expand=True)[0].str.replace('"', '')
gtf['Gene_ID']=gtf.Attributes.str.split('gene_id ', expand=True)[1].str.split(';', expand=True)[0].str.replace('"', '').str.split('.', expand=True)[0]

type_dict=dict(zip(gtf.Gene_name.to_list(), gtf.Gene_type.to_list()))

df['Gene_type']=df.Gene.apply(lambda x: ';'.join([type_dict[i] if i in type_dict.keys() else 'NF' for i in x.split(';')]))

# Drop any annotations for non-coding genes - this will also drop any genes for which we could not define a gene type
df=df[(df.Gene_type.str.contains('protein_coding'))]
print(df.Gene_type.value_counts())

# For any annotations for multiple genes, restrict to protein_coding genes
# All multi-gene annotations are "upstream/downstream"
print(df[df.Gene.str.contains(';')])
def fix_anno(row, type_dict=type_dict, gtf=gtf):
	vid=row.VID
	chrom=vid.split('_')[0]
	pos=int(vid.split('_')[1])
	genes=row.Gene.split(';')
	funcs=row.Func.split(';')

	outgenes=''
	outfuncs=''

	for gene in genes:
		if gene not in type_dict or type_dict[gene]!='protein_coding':
			continue
		gtf_line=gtf[(gtf.Gene_name==gene) & (gtf.Chr==chrom)]
		gstart=gtf_line.Start.to_list()[0]
		gend=gtf_line.End.to_list()[0]

		if gene=='GPR162':
			print(gtf_line)
			print(pos)

		for func in funcs:
			if func=='downstream':
				# Gene GPR162 annotations in GTF file do not agree with ANNOVAR annotations - default to ANNOVAR annotations
				if pos>(gend-2000) and pos<=gend+1000:
					outgenes+=f'{gene};'
					outfuncs+=f'{func};'
			if func=='upstream':
				if pos<gstart and pos>=gstart-1000:
					outgenes+=f'{gene};'
					outfuncs+=f'{func};'

	outgenes=outgenes.strip(';')
	outfuncs=outfuncs.strip(';')

	return outgenes, outfuncs

genes=[]
funcs=[]
for idx, row in df[df.Gene.str.contains(';')].iterrows():
	og, of=fix_anno(row)
	genes.append(og)
	funcs.append(of)
df.loc[df.Gene.str.contains(';'), 'Func']=funcs
df.loc[df.Gene.str.contains(';'), 'Gene']=genes

# Split any remaining genes
df.Gene=df.Gene.str.split(';')
df.Func=df.Func.str.split(';')

df=df.explode(['Gene', 'Func'])

# Annotate gene ID
id_dict=dict(zip(gtf.Gene_name.to_list(), gtf.Gene_ID.to_list()))
df['Gene_ID']=df.Gene.map(id_dict)

df=df[['Sample', 'VID', 'REF_LEN', 'ALT_LEN', 'Zygosity', 'ENSEMBL_freq', 'ENSEMBL_Mean', 'ENSEMBL_SD', 'Intracohort_count',
	   	'Func', 'Gene', 'ExonicFunc', 'Gene_ID']]

# Save
df.to_csv(f'{OUTPUT_DIR}/Full_STR_calls.csv', index=False)

# Annotate variant types
df['variant_type']='STR_'+df.Func
df['variant_type_2']='STR_noncoding'
df.loc[df.Func.isin(['exonic', 'splicing']), 'variant_type_2']='STR_coding'

df=df[['Sample', 'VID', 'Gene_ID', 'Gene', 'variant_type', 'variant_type_2']]

# Save
df.to_csv(f'{OUTPUT_DIR}/Final_STR_calls.csv', index=False)

