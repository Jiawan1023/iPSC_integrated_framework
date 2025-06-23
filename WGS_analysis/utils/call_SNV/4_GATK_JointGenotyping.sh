#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=GATK
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/WGS_processing/src # TODO: set dir to project dir
#SBATCH -o logs/4_GATK_JointGenotyping/%a.log # TODO: set slurm log file
#SBATCH -e logs/4_GATK_JointGenotyping/%a.log # TODO: set slurm log file
#SBATCH --array 1-50

# Merge GVCF by samples

echo `date` starting job on $HOSTNAME

OUT_DIR="/data7/WGS_processing/data" # TODO: set project dir path
GATK_IMAGE="/data6/deepro/computational_pipelines/gatk/gatk.sif" # TODO: set GATK pulled image path

REF_LOC="/data7/WGS_processing/src/files/hg38/" # TODO: set reference file location
REF_FILE="Homo_sapiens_assembly38.fasta"  # TODO: set reference filename

INTERVAL_FILE_LIST="files/hg38/hg38_wgs_scattered_calling_intervals.txt"
INTERVAL=`head -n $SLURM_ARRAY_TASK_ID $INTERVAL_FILE_LIST | tail -n 1`

GVCFs=`cat files/gvcf_input.txt`

# Make tmp dir
mkdir -p $OUT_DIR/tmp/JointGenotyping/$SLURM_ARRAY_TASK_ID

# Combine GVCF
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data -B $REF_LOC:/input_data $GATK_IMAGE /gatk/gatk \
CombineGVCFs \
	-V $GVCFs \
	-L /input_data/$INTERVAL \
	-R /input_data/$REF_FILE \
	-O /data/JointGenotyping/$SLURM_ARRAY_TASK_ID.g.vcf.gz \
	--tmp-dir /data/tmp/JointGenotyping/$SLURM_ARRAY_TASK_ID

# Genotype
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data -B $REF_LOC:/input_data $GATK_IMAGE /gatk/gatk \
GenotypeGVCFs \
        -R /input_data/$REF_FILE \
        -V /data/JointGenotyping/$SLURM_ARRAY_TASK_ID.g.vcf.gz \
        -L /input_data/$INTERVAL \
        -O /data/JointGenotyping/$SLURM_ARRAY_TASK_ID.scatter.g.vcf.gz \
        -G StandardAnnotation -G AS_StandardAnnotation \
        --allow-old-rms-mapping-quality-annotation-data \
        --merge-input-intervals \
        --tmp-dir /data/tmp/JointGenotyping/$SLURM_ARRAY_TASK_ID

# Clean up tmp files
rm -r $OUT_DIR/tmp/JointGenotyping/$SLURM_ARRAY_TASK_ID

echo Ended job at `date`
