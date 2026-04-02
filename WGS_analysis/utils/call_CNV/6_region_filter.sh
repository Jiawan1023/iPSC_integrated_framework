#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=gnomad_CNV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/iPSC_corrine/CNV/src
#SBATCH -o logs/2_mappability.log
#SBATCH -e logs/2_mappability.log

echo Started at `date` on $HOSTNAME

# FIlter CNV calls based on problematic regions in hg38

# Input and output
INPUT_CNV=../data/gnomAD_filter/CNV_anno/merged_cnv.bed
UCSC_DIR=files/UCSC_files
GENE_ANNO=/data6/johnathan/cnv_calling_input/data/cnv_annotation_files/use_these
OUTPUT_DIR=../data/region_filter

BEDTOOLS=/data7/software/bedtools/bedtools2/bin

mkdir $OUTPUT_DIR

# Filter SD regions
echo `date`: Filtering SD regions
mkdir -p $OUTPUT_DIR/UCSC_region_BED
mkdir -p $OUTPUT_DIR/UCSC_region_filter
gunzip -c $UCSC_DIR/genomicSuperDups.txt.gz | cut -f 2-4 > $OUTPUT_DIR/UCSC_region_BED/genomicSuperDups.bed
$BEDTOOLS/bedtools intersect -a $INPUT_CNV -b $OUTPUT_DIR/UCSC_region_BED/genomicSuperDups.bed -wa -v -f 0.5 > $OUTPUT_DIR/UCSC_region_filter/SD_filter.bed

# Filter centromeres
echo `date`: Filtering centromeres
gunzip -c $UCSC_DIR/centromeres.txt.gz | cut -f 2-4 > $OUTPUT_DIR/UCSC_region_BED/centromeres.bed
$BEDTOOLS/bedtools intersect -a $OUTPUT_DIR/UCSC_region_filter/SD_filter.bed -b $OUTPUT_DIR/UCSC_region_BED/centromeres.bed -wa -v -f 0.5 > $OUTPUT_DIR/UCSC_region_filter/centromere_filter.bed

# Filter problematic regions
echo `date`: Filtering problematic regions
gunzip -c $UCSC_DIR/UCSC_Unusual_regions_hg38.bed.gz | cut -f 1-3 > $OUTPUT_DIR/UCSC_region_BED/Unusual_regions.bed
$BEDTOOLS/bedtools intersect -a $OUTPUT_DIR/UCSC_region_filter/centromere_filter.bed -b $OUTPUT_DIR/UCSC_region_BED/Unusual_regions.bed -wa -v -f 0.5 > $OUTPUT_DIR/UCSC_region_filter/Unusual_region_filter.bed

gunzip -c $UCSC_DIR/ENCODE_Blacklist_hg38.bed.gz | cut -f 1-3 > $OUTPUT_DIR/UCSC_region_BED/ENCODE_Blacklist.bed
$BEDTOOLS/bedtools intersect -a $OUTPUT_DIR/UCSC_region_filter/Unusual_region_filter.bed -b $OUTPUT_DIR/UCSC_region_BED/ENCODE_Blacklist.bed -wa -v -f 0.5 > $OUTPUT_DIR/UCSC_region_filter/ENCODE_Blacklist_filter.bed

gunzip -c $UCSC_DIR/GRC_Exclusions_hg38.bed.gz | cut -f 1-3 > $OUTPUT_DIR/UCSC_region_BED/GRC_Exclusions.bed
$BEDTOOLS/bedtools intersect -a $OUTPUT_DIR/UCSC_region_filter/ENCODE_Blacklist_filter.bed -b $OUTPUT_DIR/UCSC_region_BED/GRC_Exclusions.bed -wa -v -f 0.5 > $OUTPUT_DIR/UCSC_region_filter/GRC_Exclusions_filter.bed

echo Ended at `date`
