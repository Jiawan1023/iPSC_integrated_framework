import pandas as pd
import scipy.stats as stats
from scipy import stats
import pybedtools

deg_filename="/Users/jks6575/Dropbox/Jiawan/Analysis/outlier_expression/Fisher_exact_atac/ipsc_A38_CC_dp_results_new.csv"
output_filename="/Users/jks6575/Dropbox/Jiawan/Analysis/outlier_expression/Fisher_exact_atac/ipsc_wgs_atac.csv"
samples=['CR001_A38', 'CR001_WT']

# Load in variant data
vardf=pd.read_csv('/Users/jks6575/Dropbox/Jiawan/Analysis/outlier_expression/variant_data/ipsc_family_variants.csv')
vardf = vardf[vardf['Sample'].isin(samples)]
# filter variant in CR001_A38 the same with CR001_WT
vardf = vardf[
    (vardf['Sample'] == 'CR001_A38') & 
    (vardf['vid'].isin(vardf[vardf['Sample'] == 'CR001_WT']['vid']))
]
print(vardf)
variant_types=[ 'upstream','downstream', 'utr5','utr3', 'intron', 
                'deletion_upstream', 'deletion_downstream', 'deletion_utr5', 'deletion_utr3', 'deletion_intron',
                'duplication_upstream', 'duplication_downstream', 'duplication_utr5', 'duplication_utr3', 'duplication_intron', 
               'STR_upstream', 'STR_downstream', 'STR_UTR5' ,'STR_UTR3', 'STR_intronic']

vardf[['chrom', 'Start', 'End']] = vardf['vid'].str.split('_', expand=True)[[0, 1, 2]]

# For cases where 'End' is empty or missing, fill 'End' with the value of 'Start'
vardf['End'] = vardf.apply(lambda row: row['Start'] if not row['End'].isdigit() else row['End'], axis=1)

print(vardf)

#variant_types=['SNV_coding', 'SNV_noncoding', 'deletion_coding', 'deletion_noncoding', 'duplication_coding', 
               #'duplication_noncoding', 'STR_coding', 'STR_noncoding']
# Load in diffpeak list
diffpeaks=pd.read_csv(deg_filename)
print(diffpeaks)

#input population peaks as peak space
columns = ['chrom', 'Start', 'End']
stat_lst=[]
peak_space=pd.read_csv('/Users/jks6575/Dropbox/Jiawan/Analysis/ATAC_seq_analysis/iPSC/ipsc_atac_peaks.bed', sep='\t', header=None, names=columns)
peak_space['peak_name'] = peak_space['chrom'] + '_' + peak_space['Start'].astype(str) + '_' + peak_space['End'].astype(str)

print(peak_space)

peak_space['Start'] = pd.to_numeric(peak_space['Start'], errors='coerce')
peak_space['End'] = pd.to_numeric(peak_space['End'], errors='coerce')

vardf['Start'] = pd.to_numeric(vardf['Start'], errors='coerce')
vardf['End'] = pd.to_numeric(vardf['End'], errors='coerce')

diffpeaks_bed = pybedtools.BedTool.from_dataframe(diffpeaks[['chrom', 'Start', 'End']])
vardf_bed = pybedtools.BedTool.from_dataframe(vardf[['chrom', 'Start', 'End','gene','gene_id','vid','variant_type']])
   
overlaps = diffpeaks_bed.intersect(vardf_bed, wa=True, wb=True)
overlap_df = overlaps.to_dataframe(names=['peak_chrom', 'peak_Start', 'peak_End', 'vardf_chrom', 'vardf_Start', 'vardf_End','gene','gene_id','vid','variant_type'])    
print(overlap_df)
overlap_df.to_csv(output_filename, index=False)



