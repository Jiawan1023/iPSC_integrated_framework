import pandas as pd

samples = ['CR001', 'CR001A38', 'CR007', 'SG037', 'SG038', 'SG039', 'SG040', 'SG041', 'SG042', 'SG043', 'SG044', 'SG045', 'SG046', 'SG213', 'SG214', 'SG215', 'SG216', 'SG217', 'SG218', 'SG219', 'SG220', 'SG221', 'SG225', 'SG226', 'SG227', 'SG378', 'SG552']
statstr_by_cohort_df = pd.read_csv('/data7/ipsc/16p12_1_del/str_calls/data/statstr/statstr_by_cohort_9-30-24.tab', sep='\t')
thousand_genome_statstr_dir = '/statstr'

final_df = pd.DataFrame()

for chrm in range(1, 23):
    print(chrm)
    thousand_genome_statstr_df = pd.read_csv(f'{thousand_genome_statstr_dir}/ensemble_chr{chrm}_filtered.tab', sep='\t')
    
    # Merge cohort and 1000G dfs
    merged_statstr_dfs = statstr_by_cohort_df.merge(thousand_genome_statstr_df.drop_duplicates(), on=['chrom','start','end'], 
                   how='left', indicator=True)

    # Add all "left_only" entries to final df
    cohort_df_only = merged_statstr_dfs[merged_statstr_dfs['_merge'] == 'left_only']
    final_df = pd.concat([final_df, cohort_df_only])

    both_dfs = merged_statstr_dfs[merged_statstr_dfs['_merge'] == 'both']

    # Expand the afreq cols
    both_dfs['afreq_x'] = both_dfs['afreq_x'].str.split(',')
    both_dfs['afreq_y'] = both_dfs['afreq_y'].str.split(',')
    both_dfs = both_dfs.explode('afreq_x')
    both_dfs = both_dfs.explode('afreq_y')
    both_dfs[['alleles_x', 'freqs_x']] = both_dfs['afreq_x'].str.split(':', expand=True)
    both_dfs[['alleles_y', 'freqs_y']] = both_dfs['afreq_y'].str.split(':', expand=True)

    # Filter for matching alleles and freq vals less than 1%
    both_dfs = both_dfs[both_dfs['alleles_x'] == both_dfs['alleles_y']]
    sig_freq_vals_df = both_dfs[both_dfs['freqs_y'].astype('float') < 0.001].reset_index(drop=True) # NOTE (2024-9-24): this was changed from 0.01 to 0.001
    final_df = pd.concat([final_df, cohort_df_only])

    non_sig_freq_vals_df = both_dfs[both_dfs['freqs_y'].astype('float') >= 0.001].reset_index(drop=True) # NOTE (2024-9-24): this was changed from 0.01 to 0.001
    non_sig_freq_vals_df = non_sig_freq_vals_df[non_sig_freq_vals_df['mean_x'] > non_sig_freq_vals_df['mean_y'] + (2 * non_sig_freq_vals_df['var_y'])]

    final_df = pd.concat([final_df, non_sig_freq_vals_df])

final_df = final_df.drop(['afreq_y', 'mean_y', 'var_y', '_merge', 'alleles_y', 'freqs_y', 'alleles_x', 'freqs_x'], axis=1)
final_df = final_df.rename(columns={'afreq_x':'afreq', 'mean_x':'mean', 'var_x':'var'})

# Add samples to statstr df
mergestr_vcf = pd.read_csv('/data7/ipsc/16p12_1_del/str_calls/data/mergestr_by_cohort_9-30-24.vcf', skiprows=range(3414), sep='\t')
mergestr_vcf = mergestr_vcf.rename(columns={'#CHROM':'chrom', 'POS':'start'})
mergestr_statstr_merge_df = final_df.merge(mergestr_vcf, on=['chrom', 'start'], how='left')

mergestr_statstr_merge_df['sample'] = ''

for col in mergestr_statstr_merge_df.columns[13:39]:
    mergestr_statstr_merge_df.loc[(mergestr_statstr_merge_df[col].apply(lambda x: x.split(':')[0]) != '0/0') & \
                                  (mergestr_statstr_merge_df[col].apply(lambda x: x.split(':')[0]) != '.'), 'sample'] += col + ','

mergestr_statstr_merge_df = mergestr_statstr_merge_df.drop(mergestr_statstr_merge_df.columns[6:39], axis=1)
final_df = mergestr_statstr_merge_df

final_df['sample'] = final_df['sample'].apply(lambda x: x.split(',')[:-1])
final_df = final_df[final_df['sample'].apply(lambda x: len(x) != 0)] # Keep rows that don't have empty sample list
# final_df = final_df.explode('sample')
print(final_df)

final_df.to_csv('/data7/ipsc/16p12_1_del/str_calls/data/filter_statstr_2024-9-30.tsv', sep='\t', index=False)
