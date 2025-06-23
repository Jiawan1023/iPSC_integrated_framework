#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=annotate_using_vep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data7/ipsc/16p12_1_del/str_calls # TODO: set dir to project dir
#SBATCH -o /data7/ipsc/16p12_1_del/str_calls/slurm/logs/6_out_%a.log # TODO: set slurm output file
#SBATCH -e /data7/ipsc/16p12_1_del/str_calls/slurm/logs/6_err_%a.log # TODO: set slurm input file
#SBATCH --exclude=sarah
#SBATCH --array 1-26

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data7/ipsc/16p12_1_del/str_calls/slurm/files/5_smap.txt)
sample_id=`echo $LINE | cut -d" " -f1`

python3 /data7/ipsc/16p12_1_del/str_calls/src/6_filter_for_str_expansion.py --sample $sample_id
