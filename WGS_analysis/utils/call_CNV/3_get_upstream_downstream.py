import pandas as pd

df = pd.read_csv('CR001.bed', sep='\t', names=['chrm', 'start', 'end', 'size', 'CNV_type', 'read_depth', 'e_val1',
        'e_val2', 'e_val3', 'e_val4', 'q0', 'pN', 'dG', 'gene_ids', 'gene_names', 
        'Func.refGene', 'info_chrm', 'info_start', 'info_end', 'info_gene_id', 'info_INSERT_WHAT_IT_IS_HERE', 'info_strand'])

upstream_df = pd.DataFrame(columns=['chrm', 'start', 'end', 'gene_id', 'INSERT_WHAT_IT_IS_HERE', 'strand'])
downstream_df = pd.DataFrame(columns=['chrm', 'start', 'end', 'gene_id', 'INSERT_WHAT_IT_IS_HERE', 'strand'])

utr3_df = pd.read_csv('utr3.bed', sep='\t', header=None)
utr5_df = pd.read_csv('utr5.bed', sep='\t', header=None)

df_size = len(df) # Since the length of the df can change in the loop, it's best to record the original length of the df here

downstream = False
for i in range(0, len(utr3_df)):
    new_row = {}
    
    if utr3_df.loc[i][5] == '+':
        range_start_dis = 1
        range_end_dis = 5000
        downstream = True
    elif utr3_df.loc[i][5] == '-':
        range_start_dis = -5000
        range_end_dis = -1
        downstream = False

    new_row = {'chrm': utr3_df.loc[i][0], 
               'start': int(utr3_df.loc[i][2]) + range_start_dis, 
               'end': int(utr3_df.loc[i][2]) + range_end_dis, 
               'gene_id': utr3_df.loc[i][3], 
               'INSERT_WHAT_IT_IS_HERE': utr3_df.loc[i][4], 
               'strand': utr3_df.loc[i][5]
               }
    
    if downstream: downstream_df.loc[len(downstream_df)] = new_row
    else: upstream_df.loc[len(upstream_df)] = new_row

for i in range(0, len(utr5_df)):
    new_row = {}
    
    if utr5_df.loc[i][5] == '-':
        range_start_dis = 1
        range_end_dis = 5000
        downstream = True
    elif utr5_df.loc[i][5] == '+':
        range_start_dis = -5000
        range_end_dis = -1
        downstream = False

    new_row = {'chrm': utr5_df.loc[i][0], 
               'start': int(utr5_df.loc[i][1]) + range_start_dis, 
               'end': int(utr5_df.loc[i][1]) + range_end_dis, 
               'gene_id': utr5_df.loc[i][3], 
               'INSERT_WHAT_IT_IS_HERE': utr5_df.loc[i][4], 
               'strand': utr5_df.loc[i][5]
               }
    
    if downstream: downstream_df.loc[len(downstream_df)] = new_row
    else: upstream_df.loc[len(upstream_df)] = new_row

upstream_df.to_csv('upstream.bed', sep='\t', index=False, header=False)
downstream_df.to_csv('downstream.bed', sep='\t', index=False, header=False)
