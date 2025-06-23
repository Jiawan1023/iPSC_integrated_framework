#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=GATK
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=4
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=30G
#SBATCH --chdir /data7/WGS_processing/src # TODO: set dir to project dir
#SBATCH -o logs/3_GATK_HaplotypeCaller/%a.log # TODO: set slurm log file
#SBATCH -e logs/3_GATK_HaplotypeCaller/%a.log # TODO: set slurm log file
#SBATCH --array 5-6
#SBATCH --exclude qingyu

# Merge GVCF by samples

echo `date` starting job on $HOSTNAME

OUT_DIR="/data7/WGS_processing/data" # TODO: set project dir path
GATK_IMAGE="/data6/deepro/computational_pipelines/gatk/gatk.sif" # TODO: set GATK pulled image path

REF_LOC="/data7/WGS_processing/src/files/hg38/" # TODO: set reference file location
REF_FILE="Homo_sapiens_assembly38.fasta"  # TODO: set reference filename

INTERVAL_FILE_LIST="files/hg38/hg38_wgs_scattered_calling_intervals.txt"

LINE=`head -n $SLURM_ARRAY_TASK_ID files/sample_info.csv | tail -n 1`

SAMPLE=`echo $LINE | cut -f 1 -d ,`

echo $SAMPLE

# Call SNVs and INDELs in parallel
mkdir -p $OUT_DIR/$SAMPLE/gvcf
for i in {1..50};do
	bash files/scripts/GATK_HaplotypeCaller.sh $SAMPLE $i $OUT_DIR $GATK_IMAGE $REF_LOC $REF_FILE $INTERVAL_FILE_LIST &
done
wait

# Merge VCFs into single file
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data $GATK_IMAGE /gatk/gatk \
MergeVcfs \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.1.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.2.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.3.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.4.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.5.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.6.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.7.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.8.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.9.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.10.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.11.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.12.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.13.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.14.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.15.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.16.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.17.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.18.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.19.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.20.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.21.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.22.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.23.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.24.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.25.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.26.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.27.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.28.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.29.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.30.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.31.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.32.g.vcf \
	--INPUT	/data/$SAMPLE/gvcf/$SAMPLE.33.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.34.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.35.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.36.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.37.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.38.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.39.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.40.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.41.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.42.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.43.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.44.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.45.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.46.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.47.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.48.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.49.g.vcf \
	--INPUT /data/$SAMPLE/gvcf/$SAMPLE.50.g.vcf \
	--OUTPUT /data/$SAMPLE/$SAMPLE.g.vcf \
	--TMP_DIR /data/tmp/$SAMPLE

# Bgzip and index
bgzip /data7/WGS_processing/data/$SAMPLE/$SAMPLE.g.vcf
tabix -p vcf /data7/WGS_processing/data/$SAMPLE/$SAMPLE.g.vcf.gz

echo Ended job at `date`
