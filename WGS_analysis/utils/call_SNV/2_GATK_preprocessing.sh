#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=GATK
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data7/WGS_processing/src # TODO: set dir to project dir
#SBATCH -o logs/2_GATK_preprocessing/%a.log # TODO: set slurm log file
#SBATCH -e logs/2_GATK_preprocessing/%a.log # TODO: set slurm log file
#SBATCH --array 5-8

# Align FASTQ using GATK Best Practices

echo `date` starting job on $HOSTNAME

OUT_DIR="/data7/WGS_processing/data" # TODO: set project dir path
GATK_IMAGE="/data6/deepro/computational_pipelines/gatk/gatk.sif" # TODO: set GATK pulled image path

REF_LOC="/data7/WGS_processing/src/files/hg38/" # TODO: set reference file location
REF_FILE="Homo_sapiens_assembly38.fasta"  # TODO: set reference filename

DBSNP_VCF="Homo_sapiens_assembly38.dbsnp138.vcf"
MG_VCF="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
INDELS_VCF="Homo_sapiens_assembly38.known_indels.vcf.gz"

INTERVAL_FILE="files/sequence_groupings/sequence_grouping.txt"

LINE=`head -n $SLURM_ARRAY_TASK_ID files/sample_info.csv | tail -n 1`

SAMPLE=`echo $LINE | cut -f 1 -d ,`
FASTQ_LOC=`echo $LINE | cut -f 2 -d ,`
R1=`echo $LINE | cut -f 3 -d ,`
R2=`echo $LINE | cut -f 4 -d ,`
RG=`echo $LINE | cut -f 5 -d ,`
LN=`echo $LINE | cut -f 6 -d ,`
PU=`echo $LINE | cut -f 7 -d ,`
PN=`echo $LINE | cut -f 8 -d ,`
SC=`echo $LINE | cut -f 9 -d ,`

echo $SAMPLE

mkdir -p /data7/WGS_processing/data/$SAMPLE
mkdir -p /data7/WGS_processing/data/tmp/$SAMPLE

# FASTQ to UBAM
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data -B $FASTQ_LOC:/input_data $GATK_IMAGE /gatk/gatk \
FastqToSam \
	--FASTQ /input_data/$R1 \
	--FASTQ2 /input_data/$R2 \
	--OUTPUT /data/$SAMPLE/$SAMPLE.all.ubam.bam \
	--READ_GROUP_NAME $RG \
	--SAMPLE_NAME $SAMPLE \
	--LIBRARY_NAME $LN \
	--PLATFORM_UNIT $PU \
	--PLATFORM $PN \
	--SEQUENCING_CENTER $SC \
	--TMP_DIR /data/tmp/$SAMPLE

# Index new bam file
samtools index $OUT_DIR/$SAMPLE/$SAMPLE.all.ubam.bam

# Mark adapter sequences
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data $GATK_IMAGE /gatk/gatk \
MarkIlluminaAdapters \
	I=/data/$SAMPLE/$SAMPLE.all.ubam.bam \
	O=/data/$SAMPLE/$SAMPLE.markadapters.bam \
	M=/data/$SAMPLE/$SAMPLE.markadapters_metrics.txt \
	TMP_DIR=/data/tmp/$SAMPLE

# Convert back to FASTQ
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data $GATK_IMAGE /gatk/gatk \
SamToFastq \
	INPUT=/data/$SAMPLE/$SAMPLE.markadapters.bam \
	FASTQ=/data/$SAMPLE/$SAMPLE.sam2fastq.fastq \
	INTERLEAVE=true \
	NON_PF=true \
	TMP_DIR=/data/tmp/$SAMPLE

# Align with BWA
/data7/software/bwa/bwa-0.7.17/bwa mem -M -t 8 -p $REF_LOC$REF_FILE $OUT_DIR/$SAMPLE/$SAMPLE.sam2fastq.fastq | samtools view -b -o $OUT_DIR/$SAMPLE/$SAMPLE.aligned.bam

# Apply and adjust metadata
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data -B $REF_LOC:/input_data $GATK_IMAGE /gatk/gatk \
MergeBamAlignment \
	--VALIDATION_STRINGENCY SILENT \
	--EXPECTED_ORIENTATIONS FR \
	--ATTRIBUTES_TO_RETAIN X0 \
	--ALIGNED_BAM /data/$SAMPLE/$SAMPLE.aligned.bam \
	--UNMAPPED_BAM /data/$SAMPLE/$SAMPLE.all.ubam.bam \
	--OUTPUT /data/$SAMPLE/$SAMPLE.mba.bam \
	--REFERENCE_SEQUENCE /input_data/$REF_FILE \
	--PAIRED_RUN true \
	--SORT_ORDER "unsorted" \
	--IS_BISULFITE_SEQUENCE false \
	--ALIGNED_READS_ONLY false \
	--CLIP_ADAPTERS false \
	--MAX_RECORDS_IN_RAM 2000000 \
	--ADD_MATE_CIGAR true \
	--MAX_INSERTIONS_OR_DELETIONS -1 \
	--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
	--PROGRAM_RECORD_ID "bwamem" \
	--PROGRAM_GROUP_VERSION "0.7.17" \
	--PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -v 3 -t 16 -Y $REF" \
	--PROGRAM_GROUP_NAME "bwamem" \
	--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
	--ALIGNER_PROPER_PAIR_FLAGS true \
	--UNMAP_CONTAMINANT_READS true \
	--TMP_DIR /data/tmp/$SAMPLE

# Mark Duplicates
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data $GATK_IMAGE /gatk/gatk \
MarkDuplicates \
	--INPUT /data/$SAMPLE/$SAMPLE.mba.bam \
	--OUTPUT /data/$SAMPLE/$SAMPLE.markdup.bam \
	--METRICS_FILE /data/$SAMPLE/$SAMPLE.duplicate_metrics.txt \
	--VALIDATION_STRINGENCY SILENT \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
	--ASSUME_SORT_ORDER "queryname" \
	--CREATE_MD5_FILE true \
	--TMP_DIR /data/tmp/$SAMPLE

# Sort and Fix Tags
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data -B $REF_LOC:/input_data $GATK_IMAGE /gatk/gatk \
SortSam \
	--INPUT /data/$SAMPLE/$SAMPLE.markdup.bam \
	--OUTPUT /data/$SAMPLE/$SAMPLE.sorted.bam \
	--SORT_ORDER "coordinate" \
	--CREATE_INDEX false \
	--CREATE_MD5_FILE false \
	--TMP_DIR /data/tmp/$SAMPLE

singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data -B $REF_LOC:/input_data $GATK_IMAGE /gatk/gatk \
SetNmMdAndUqTags \
	--INPUT /data/$SAMPLE/$SAMPLE.sorted.bam \
	--OUTPUT /data/$SAMPLE/$SAMPLE.sortfix.bam \
	--CREATE_INDEX true \
	--CREATE_MD5_FILE true \
	--REFERENCE_SEQUENCE /input_data/$REF_FILE \
	--TMP_DIR /data/tmp/$SAMPLE

# Calculate Base Quality Score Recalibration in parallel
mkdir -p $OUT_DIR/$SAMPLE/recalibration_reports
for i in {1..18}; do
	bash files/scripts/GATK_BQSR.sh $SAMPLE $i $OUT_DIR $GATK_IMAGE $REF_LOC $REF_FILE $DBSNP_VCF $MG_VCF $INDELS_VCF $INTERVAL_FILE &
done
wait

# Compile BQSR reports
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data $GATK_IMAGE /gatk/gatk \
GatherBQSRReports \
	-I /data/$SAMPLE/recalibration_reports/1.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/2.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/3.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/4.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/5.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/6.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/7.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/8.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/9.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/10.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/11.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/12.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/13.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/14.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/15.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/16.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/17.recal_data.csv \
	-I /data/$SAMPLE/recalibration_reports/18.recal_data.csv \
	-O /data/$SAMPLE/$SAMPLE.recal_data.csv \
	--tmp-dir /data/tmp/$SAMPLE

# Apply BQSR in parallel
for i in {1..18}; do
	bash files/scripts/GATK_ApplyBQSR.sh $SAMPLE $i $OUT_DIR $GATK_IMAGE $REF_LOC $REF_FILE $INTERVAL_FILE &
done
wait

# Gather BAMs from BQSR
singularity exec --containall -H $OUT_DIR:/data -B $OUT_DIR:/data $GATK_IMAGE /gatk/gatk \
GatherBamFiles \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.1.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.2.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.3.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.4.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.5.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.6.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.7.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.8.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.9.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.10.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.11.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.12.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.13.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.14.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.15.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.16.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.17.bam \
	--INPUT /data/$SAMPLE/$SAMPLE.BQSR.18.bam \
	--OUTPUT /data/$SAMPLE/$SAMPLE.bam \
	--CREATE_INDEX true \
	--CREATE_MD5_FILE true \
	--TMP_DIR /data/tmp/$SAMPLE

# Index new bam file
samtools index $OUT_DIR/$SAMPLE/$SAMPLE.bam

echo Ended job at `date`
