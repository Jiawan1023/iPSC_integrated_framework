import pandas as pd
import argparse

def main(file):
    df = pd.read_csv(file, sep='\t')
    df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT']]
    df.to_csv(file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='4_filter_data')
    parser.add_argument("--file", type=str)

    cli_args = parser.parse_args()

    main(cli_args.chr)
