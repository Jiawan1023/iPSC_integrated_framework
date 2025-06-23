import os
import cnvpytor
import argparse
import pandas as pd


def call_cnvs(sample_name, save_dir, binsize, max_threads=64):
    # create pytor file and instance
    pytor_file = os.path.join(save_dir, "pytors", f"{sample_name}.pytor")
    app = cnvpytor.Root(pytor_file, create=False, max_cores=max_threads)
    app.calculate_histograms([binsize])
    app.partition([binsize])
    calls = app.call([binsize])
    colnames = ["CNV_type", "chrm", "start", "end", "size", "read_depth", "e_val1", "e_val2", "e_val3", "e_val4", "q0", "pN", "dG"]
    df = pd.DataFrame(data=calls[binsize], columns=colnames)
    call_file = os.path.join(save_dir, "calls", f"{sample_name}_{binsize}.tsv")
    df.to_csv(call_file, index=False, sep="\t")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CNVpytor caller")
    # required arguments
    parser.add_argument("--sample_name", type=str, help="name of the sample, all files will have this as prefix")
    parser.add_argument("--save_dir", type=str, help="Path where results will be saved")
    parser.add_argument("--bin_size", type=str, help="bin size")
    cli_args = parser.parse_args()
    call_cnvs(cli_args.sample_name, cli_args.save_dir, int(cli_args.bin_size))
