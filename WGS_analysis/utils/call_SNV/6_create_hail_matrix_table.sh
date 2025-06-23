#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=HAIL_matrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/WGS_processing/src # TODO: set dir to cache dir
#SBATCH -o logs/6_create_matrix.log # TODO: set slurm output file
#SBATCH -e logs/6_create_matrix.log # TODO: set slurm input file
#SBATCH --nodelist=sarah

echo `date` starting job on $HOSTNAME

export HOME="/data7/WGS_processing"
source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/hail_env

OUTDIR=../data/Isogenic_lines/HAIL

mkdir -p $OUTDIR

python files/scripts/HAIL_create_matrix.py ../data/Isogenic_lines/allsample.allchr.merged.g.vcf.gz $OUTDIR/allsample.hail.mt GRCh38

# Move the logfile it makes out of the main directory
mv hail*.log $OUTDIR/allsample.hail.mt

echo `date` ending job on $HOSTNAME
