import os
import cnvpytor
import logging
import argparse


def create_pytor(cram_file, sample_name, save_dir, max_threads=64):
    # define log file
    logfile = os.path.join(save_dir, "logs", f"{sample_name}.log")
    print(logfile)
    logging.basicConfig(
        filename=logfile, level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
    logger = logging.getLogger('cnvpytor')
    # create pytor file and instance
    pytor_file = os.path.join(save_dir, "pytors", f"{sample_name}.pytor")
    app = cnvpytor.Root(pytor_file, create=True, max_cores=max_threads)
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    # import and save rd signal
    app.rd([cram_file], chroms=chroms)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CNVpytor caller")
    # required arguments
    parser.add_argument("--cram_file", type=str, help="path to the cram file")
    parser.add_argument("--sample_name", type=str, help="name of the sample, all files will have this as prefix")
    parser.add_argument("--save_dir", type=str, help="Path where results will be saved")
    
    cli_args = parser.parse_args()
    create_pytor(cli_args.cram_file, cli_args.sample_name, cli_args.save_dir)
