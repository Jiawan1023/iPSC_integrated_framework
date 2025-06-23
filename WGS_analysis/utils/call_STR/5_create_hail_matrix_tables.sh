#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=create_hail_matrix_tables
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data7/ipsc/16p12_1_del/str_calls # TODO: set dir to project dir
#SBATCH -o /data7/ipsc/16p12_1_del/str_calls/slurm/logs/5_out_%a.log # TODO: set slurm output file
#SBATCH -e /data7/ipsc/16p12_1_del/str_calls/slurm/logs/5_err_%a.log # TODO: set slurm input file
#SBATCH --nodelist=sarah

export HOME="/data7/ipsc/16p12_1_del/str_calls"

source /opt/anaconda/bin/activate /data6/johnathan/miniconda3/miniconda/envs/hail4

python3 /data7/ipsc/16p12_1_del/str_calls/src/5_create_hail_matrix_tables.py
