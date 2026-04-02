#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=iPSC_STR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/iPSC_corrine/STR/src
#SBATCH -o logs/1_dumpSTR.log
#SBATCH -e logs/1_dumpSTR.log

echo Started at `date` on $HOSTNAME

export HOME="/data7/iPSC_corrine"
source /opt/anaconda/bin/activate /data6/deepro/miniconda3/envs/strenv

# Run dumpSTR to remove low quality sites and call CNVs on filtered sites

# Input and output
MERGESTR=/data7/ipsc/16p12_1_del/str_calls/data/mergestr_by_cohort/mergestr_by_cohort_9-30-24.vcf
SUPDUP=files/genomicSuperDups_hg38.bed # Created from UCSC Table Browser as in WGS paper
ENSEMBLTR=/data7/johnathan/ensembleTR/statstr
SMAP='/data7/ipsc/16p12_1_del/str_calls/slurm/files/1_smap.txt'
REF='/data7/WGS_processing/src/files/hg38/Homo_sapiens_assembly38.fasta'
OUTPUT_DIR=../data/Locus_filter

BCFTOOLS=/data5/anastasia/sw/bcftools-1.12

mkdir -p $OUTPUT_DIR

# Compress input VCF
echo `date`: Compressing VCF
bgzip -c $MERGESTR > $OUTPUT_DIR/mergestr_by_cohort_9-30-24.vcf.gz
tabix -p vcf -f $OUTPUT_DIR/mergestr_by_cohort_9-30-24.vcf.gz

# Filter SuperDup file
echo `date`: Filtering SuperDups
cut -f 2-5 $SUPDUP | tail -n+2 | sort -k1,1V -k2,2n | bgzip > files/genomicSuperDups_hg38_sorted.bed.bgz
tabix -p bed -f files/genomicSuperDups_hg38_sorted.bed.bgz

# Run dumpSTR
echo `date`: Running dumpSTR
dumpSTR \
	--vcf $OUTPUT_DIR/mergestr_by_cohort_9-30-24.vcf.gz \
	--out $OUTPUT_DIR/1_dumpSTR \
	--min-locus-hwep 0.00001 \
	--min-locus-callrate 0.8 \
	--filter-regions files/genomicSuperDups_hg38_sorted.bed.bgz  \
	--filter-regions-names SEGDUP

bgzip $OUTPUT_DIR/1_dumpSTR.vcf
tabix -p vcf -f $OUTPUT_DIR/1_dumpSTR.vcf.gz

# Extract sites passing QC
echo `date`: Extracting sites
$BCFTOOLS/bcftools view -i 'FILTER="PASS"' $OUTPUT_DIR/1_dumpSTR.vcf.gz | $BCFTOOLS/bcftools query -f '%CHROM\t%POS\t%INFO/END\n' > $OUTPUT_DIR/pass_sites.bed

# Extract loci present in ensemblTR
echo `date`: Extracting autosomal Ensembl sites
rm -rf $OUTPUT_DIR/ensemble_sites.bed
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	cut -f 1-3 $ENSEMBLTR/ensemble_chr${chr}_filtered.tab | tail -n+2 >> $OUTPUT_DIR/ensemble_sites.bed
done

# Remove QC-failed loci from calls, exclude loci not present in Ensembl, extract non-ref loci, and convert to bed
echo `date`: Removing sites
mkdir -p $OUTPUT_DIR/filter_loci
$BCFTOOLS/bcftools view -R $OUTPUT_DIR/pass_sites.bed -i 'GT[*]="alt"' $OUTPUT_DIR/mergestr_by_cohort_9-30-24.vcf.gz | bgzip > $OUTPUT_DIR/filter_loci/alt_alleles.vcf.gz
tabix -p vcf -f $OUTPUT_DIR/filter_loci/alt_alleles.vcf.gz
$BCFTOOLS/bcftools view -R $OUTPUT_DIR/ensemble_sites.bed $OUTPUT_DIR/filter_loci/alt_alleles.vcf.gz | bgzip > $OUTPUT_DIR/filter_loci/ensembl_sites.vcf.gz
tabix -p vcf -f $OUTPUT_DIR/filter_loci/ensembl_sites.vcf.gz
$BCFTOOLS/bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/PERIOD\t%INFO/RU\n' $OUTPUT_DIR/filter_loci/ensembl_sites.vcf.gz > $OUTPUT_DIR/STR_QC.bed

echo Ended at `date`
