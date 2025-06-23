import pandas as pd
import argparse


def main(sample):
    # df = pd.read_csv(f'/data7/ipsc/16p12_1_del/str_calls/data/filter_with_annovar/{sample}/{sample}.hg38_multianno.vcf', sep='\t')
    # start =  df.loc[df.FILE-START == '#CHROM'].index[0]
    # df = pd.read_csv(f'/data7/ipsc/16p12_1_del/str_calls/data/filter_with_annovar/{sample}/{sample}.hg38_multianno.vcf', sep='\t', skiprows = start)
    # print(df)
    # return

    with open(f'/data7/ipsc/16p12_1_del/str_calls/data/filter_with_annovar/{sample}/{sample}.hg38_multianno.vcf', 'r') as fh:
        rows_to_skip = 0
        while line := fh.readline():
            if line.find('#CHROM') != -1:
                break
            rows_to_skip += 1
    
    df = pd.read_csv(f'/data7/ipsc/16p12_1_del/str_calls/data/filter_with_annovar/{sample}/{sample}.hg38_multianno.vcf', sep='\t', skiprows=rows_to_skip)
    
    return_df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']]
    
    func_refgene_list = []
    gene_refgene_list = []
    for i in range(0, len(df)):
        info_list = str(df['INFO'][i]).split(';')

        func_refgene_list.append(info_list[6].split('=')[1].replace('\\x3b', ';'))
        gene_refgene_list.append(info_list[7].split('=')[1].replace('\\x3b', ';'))

    return_df.insert(len(return_df.columns), 'func_refgene', func_refgene_list)
    return_df.insert(len(return_df.columns), 'gene_refgene', gene_refgene_list)
    
    return return_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read VCF')
    parser.add_argument("--sample", type=str)

    cli_args = parser.parse_args()
    
    df = main(cli_args.sample)
    print(df)
    df.to_csv(f'/data7/ipsc/16p12_1_del/str_calls/data/filter_with_annovar/{cli_args.sample}.vcf', sep='\t')
