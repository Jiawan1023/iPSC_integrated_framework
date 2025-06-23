#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=HAIL_table
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=200G
#SBATCH --chdir /data7/johnathan/
#SBATCH -o /data7/johnathan/test_logs/%a.log # TODO: set slurm output file
#SBATCH -e /data7/johnathan/test_logs/%a.log # TODO: set slurm input file
#SBATCH --array 1-24
#SBATCH --nodelist sarah

echo `date` starting job on $HOSTNAME

REF="GRCh38"
NEW_TMP="/data7/WGS_processing/data/tmp"
INDIR="/data7/WGS_processing/data/HAIL/VEP_anno/VEP"

export HOME="/data7/WGS_processing"
source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/hail_env

export SPARK_LOCAL_DIRS=$NEW_TMP
export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=/data7/WGS_processing/data/tmp -Djava.util.prefs.systemRoot=/data7/WGS_processing/data/tmp -Duser.home=/data7/WGS_processing/data/tmp"

chr=`head -n $SLURM_ARRAY_TASK_ID /data7/WGS_processing/src/files/chromosomes.list | tail -n 1`
OUTPUT="/data7/WGS_processing/data/HAIL/tables/"$chr".csv"

python /data7/WGS_processing/src/files/scripts/HAIL_2table.py $INDIR $OUTPUT $NEW_TMP $REF $chr

echo `date` ending job on $HOSTNAME
