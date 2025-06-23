#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=VEP_annotation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=200G
#SBATCH --chdir /data7/WGS_processing/src # TODO: set dir to cache dir
#SBATCH -o logs/7_VEP_annotation/%a.log # TODO: set slurm output file
#SBATCH -e logs/7_VEP_annotation/%a.log # TODO: set slurm input file
#SBATCH --array 7

echo `date` starting job on $HOSTNAME

REF="GRCh38"
NEW_TMP="../data/tmp"
OUTDIR="../data/HAIL/VEP_anno/VEP"

export HOME="/data7/WGS_processing"
source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/hail_env

export SPARK_LOCAL_DIRS=$NEW_TMP
export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=/data7/WGS_processing/data/tmp -Djava.util.prefs.systemRoot=/data7/WGS_processing/data/tmp -Duser.home=/data7/WGS_processing/data/tmp"

chr=`head -n $SLURM_ARRAY_TASK_ID files/chromosomes.list | tail -n 1`

python files/scripts/HAIL_VEP_annotation.py ../data/HAIL/allsample.hail.mt $OUTDIR $NEW_TMP $REF $chr

echo `date` ending job on $HOSTNAME
