#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=filter_data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data7/ipsc/16p12_1_del/str_calls # TODO: set dir to project dir
#SBATCH -o /data7/ipsc/16p12_1_del/str_calls/slurm/logs/4_out_%a.log # TODO: set slurm output file
#SBATCH -e /data7/ipsc/16p12_1_del/str_calls/slurm/logs/4_err_%a.log # TODO: set slurm input file
#SBATCH --exclude=durga,ramona,laila # TODO: set nodelist
#SBATCH --array 2-26%4

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data7/ipsc/16p12_1_del/str_calls/slurm/files/4_smap.txt)
sample_id=`echo $LINE | cut -d" " -f1`
vcf_files=`echo $LINE | cut -d " " -f2`
outdir="/data7/ipsc/16p12_1_del/str_calls/data/filtered_data/${sample_id}"
mkdir -p $outdir
outpre="${outdir}/${sample_id}.vcf"

sed '1,3389d' ${vcf_files} > ${outpre}
python3 /data7/ipsc/16p12_1_del/str_calls/src/4_filter_data.py --file ${outpre}
