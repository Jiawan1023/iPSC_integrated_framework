#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=mergestr 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data7/ipsc/16p12_1_del/str_calls # TODO: set dir to project dir
#SBATCH -o /data7/ipsc/16p12_1_del/str_calls/slurm/logs/3_out_%a.log # TODO: set slurm output file
#SBATCH -e /data7/ipsc/16p12_1_del/str_calls/slurm/logs/3_err_%a.log # TODO: set slurm input file
#SBATCH --exclude=durga,ramona,laila # TODO: set nodelist
#SBATCH --array 3-26%4

export HOME="/data7/ipsc/16p12_1_del/str_calls"

source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/strenv

echo `date` starting job on $HOSTNAME

LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data7/ipsc/16p12_1_del/str_calls/slurm/files/3_smap.txt)
sample_id=`echo $LINE | cut -d" " -f1`
vcf_files=`echo $LINE | cut -d " " -f2`
refgenome="/data7/WGS_processing/src/files/hg38/Homo_sapiens_assembly38.fasta"
tr_coords="/data7/ipsc/16p12_1_del/str_calls/data/hg38_ver13.bed"
outdir="/data7/ipsc/16p12_1_del/str_calls/data/mergestr/${sample_id}"
mkdir -p $outdir
outpre="${outdir}/${sample_id}"

echo $sample_id
echo $bam_files

bgzip $vcf_files
tabix -p vcf ${vcf_files}.gz

mergeSTR \
  --vcfs ${vcf_files}.gz \
  --out $outpre
  
echo `date` ending job on $HOSTNAME
