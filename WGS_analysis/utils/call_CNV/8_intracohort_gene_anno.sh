#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=gnomad_CNV
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/iPSC_corrine/CNV/src
#SBATCH -o logs/4_intracohort_genic.log
#SBATCH -e logs/4_intracohort_genic.log

echo Started at `date` on $HOSTNAME

# FIlter CNV calls based on rarity in gnomAD

# Input and output
CNV_DIR=../data/CNV_merge
GENE_ANNO=/data6/johnathan/cnv_calling_input/data/cnv_annotation_files/use_these
OUTPUT_DIR=../data/intracohort_gene_filter

BEDTOOLS=/data7/software/bedtools/bedtools2/bin

# Annotate CNVs with intracohort frequency
echo `date`: Annotating intracohort frequency
mkdir -p $OUTPUT_DIR/Intracohort_anno
$BEDTOOLS/bedtools intersect -a $CNV_DIR/del_merged.bed -b $CNV_DIR/del_merged.bed -loj -wa -c -f 0.5 -r > $OUTPUT_DIR/Intracohort_anno/del_anno.bed
$BEDTOOLS/bedtools intersect -a $CNV_DIR/dup_merged.bed -b $CNV_DIR/dup_merged.bed -loj -wa -c -f 0.5 -r > $OUTPUT_DIR/Intracohort_anno/dup_anno.bed

# Concat DEL and DUP
echo `date`: Concatting DEL and DUP
cat $OUTPUT_DIR/Intracohort_anno/del_anno.bed > $OUTPUT_DIR/Intracohort_anno/cnv_anno.bed
cat $OUTPUT_DIR/Intracohort_anno/dup_anno.bed >> $OUTPUT_DIR/Intracohort_anno/cnv_anno.bed

# Annotate overlaps with gene features
echo `date`: Annotating features
mkdir -p $OUTPUT_DIR/Gene_anno
for feat in upstream utr5 intron exon utr3 downstream
do
	$BEDTOOLS/bedtools intersect -a $OUTPUT_DIR/Intracohort_anno/cnv_anno.bed -b $GENE_ANNO/$feat.bed -wa -wb -loj >> $OUTPUT_DIR/Gene_anno/${feat}_annotated.bed
done

echo Ended at `date`
