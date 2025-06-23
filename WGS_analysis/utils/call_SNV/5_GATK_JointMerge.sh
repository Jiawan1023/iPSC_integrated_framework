#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=GATK
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/WGS_processing/src # TODO: set dir to project dir
#SBATCH -o logs/5_GATK_JointMerge.log # TODO: set slurm log file
#SBATCH -e logs/5_GATK_JointMerge.log # TODO: set slurm log file

# Merge GVCF by samples

echo `date` starting job on $HOSTNAME

OUT_DIR="/data7/WGS_processing/data" # TODO: set project dir path
GATK_IMAGE="/data6/deepro/computational_pipelines/gatk/gatk.sif" # TODO: set GATK pulled image path

# Merge VCFs into single file
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data $GATK_IMAGE /gatk/gatk \
MergeVcfs \
	--INPUT /data/JointGenotyping/1.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/2.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/3.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/4.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/5.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/6.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/7.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/8.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/9.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/10.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/11.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/12.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/13.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/14.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/15.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/16.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/17.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/18.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/19.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/20.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/21.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/22.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/23.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/24.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/25.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/26.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/27.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/28.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/29.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/30.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/31.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/32.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/33.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/34.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/35.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/36.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/37.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/38.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/39.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/40.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/41.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/42.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/43.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/44.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/45.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/46.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/47.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/48.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/49.scatter.g.vcf.gz \
	--INPUT /data/JointGenotyping/50.scatter.g.vcf.gz \
	--OUTPUT /data/Isogenic_lines/allsample.allchr.merged.g.vcf.gz \
	--TMP_DIR /data/tmp

echo Ended job at `date`
