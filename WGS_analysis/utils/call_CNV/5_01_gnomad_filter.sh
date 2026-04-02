#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=gnomad_CNV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/iPSC_corrine/CNV/src
#SBATCH -o logs/1_gnomad.log
#SBATCH -e logs/1_gnomad.log

echo Started at `date` on $HOSTNAME

# FIlter CNV calls based on rarity in gnomAD

# Input and output
CNV_DIR=/data6/johnathan/cnv_calling_input/data/calls_with_gene_names_bed
SMAP='/data7/ipsc/16p12_1_del/str_calls/slurm/files/1_smap.txt'
GNOMAD_SITES=/data6/johnathan/chr18/remove_gnomad_cnv/data/gnomad.v4.1.sv.sites.bed
GENE_ANNO=/data6/johnathan/cnv_calling_input/data/cnv_annotation_files/use_these
OUTPUT_DIR=../data/gnomAD_filter

BEDTOOLS=/data7/software/bedtools/bedtools2/bin

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# Concat CNV calls
echo `date`: Combining CNV calls
mkdir -p Cohort_BED
rm -rf Cohort_BED/*
while read line
do
	SAMPLE=`echo $line | cut -f 1 -d ' '`
	BED=$CNV_DIR/$SAMPLE.bed

	grep "deletion" $BED | sed "s/$/\t$SAMPLE/" >> Cohort_BED/del.bed
	grep "duplication" $BED | sed "s/$/\t$SAMPLE/" >> Cohort_BED/dup.bed

done < $SMAP

# Parse gnomAD SV sites
echo `date`: Parsing gnomAD SV
mkdir -p gnomadSV_sites
grep "_DEL_" $GNOMAD_SITES | cut -f 1-4,49 | tail -n+2 | awk '{ if ($5 > 0.01) print $0}' > gnomadSV_sites/gnomadSV_common_del.bed
grep "_DUP_" $GNOMAD_SITES | cut -f 1-4,49 | tail -n+2 | awk '{ if ($5 > 0.01) print $0}' > gnomadSV_sites/gnomadSV_common_dup.bed

grep "_DEL_" $GNOMAD_SITES | cut -f 1-4,49 | tail -n+2 | awk '{ if ($5 <= 0.01) print $0}' > gnomadSV_sites/gnomadSV_rare_del.bed
grep "_DUP_" $GNOMAD_SITES | cut -f 1-4,49 | tail -n+2 | awk '{ if ($5 <= 0.01) print $0}' > gnomadSV_sites/gnomadSV_rare_dup.bed

# Filter CNV calls by gnomAD SV
echo `date`: Filtering calls
mkdir -p filtered_calls
$BEDTOOLS/bedtools intersect -a Cohort_BED/del.bed -b gnomadSV_sites/gnomadSV_common_del.bed -wa -f 0.5 -r -v > filtered_calls/del_rare.bed
$BEDTOOLS/bedtools intersect -a Cohort_BED/dup.bed -b gnomadSV_sites/gnomadSV_common_dup.bed -wa -f 0.5 -r -v > filtered_calls/dup_rare.bed

# Annotate with gnomAD frequency
echo `date`: Annotating frequency
mkdir -p CNV_anno
$BEDTOOLS/bedtools intersect -a filtered_calls/del_rare.bed -b gnomadSV_sites/gnomadSV_rare_del.bed -wa -f 0.5 -r -wb -loj | cut -f 1-16,21 > CNV_anno/gnomad_freq_del.bed
$BEDTOOLS/bedtools intersect -a filtered_calls/dup_rare.bed -b gnomadSV_sites/gnomadSV_rare_dup.bed -wa -f 0.5 -r -wb -loj | cut -f 1-16,21 > CNV_anno/gnomad_freq_dup.bed

# Merge calls
echo `date`: Merging calls
cat CNV_anno/gnomad_freq_del.bed > CNV_anno/merged_cnv.bed
cat CNV_anno/gnomad_freq_dup.bed >> CNV_anno/merged_cnv.bed

echo Ended at `date`
