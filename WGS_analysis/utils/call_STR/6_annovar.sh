#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=ANNOVAR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data7/iPSC_corrine/STR/src
#SBATCH -o logs/4_ANNOVAR.log
#SBATCH -e logs/4_ANNOVAR.log
#SBATCH --nodelist sarah

# Run ANNOVAR on STR variants to identify gene impacts

echo Started on $HOSTNAME at `date`

# Input and output
VCF_DIR=../data/GangSTR
FILTERED_STR=../data/Filtered_STRs/Filtered_loci.csv
OUTPUT_DIR=../data/ANNOVAR

BCFTOOLS=/data5/anastasia/sw/bcftools-1.12
ANNOVAR=/data4/software/annovar

mkdir -p $OUTPUT_DIR

# Merge chromosome VCFs into single file and filter for loci that passed filtering
echo `date`: Combining VCFs
cut -f 2 -d , $FILTERED_STR | tail -n+2 | sort | uniq | sed 's/:/\t/' > $OUTPUT_DIR/Filtered_sites.txt
$BCFTOOLS/bcftools concat -a -R $OUTPUT_DIR/Filtered_sites.txt $VCF_DIR/chr{1..22}/GangSTR.vcf.gz | bgzip > $OUTPUT_DIR/All_chr.GangSTR.vcf.gz
tabix -p vcf $OUTPUT_DIR/All_chr.GangSTR.vcf.gz

# Split multiallelic sites
echo `date`: Normalizing VCF
$BCFTOOLS/bcftools norm  -m-any -N $OUTPUT_DIR/All_chr.GangSTR.vcf.gz | bgzip > $OUTPUT_DIR/Normalized.vcf.gz
tabix -p vcf $OUTPUT_DIR/Normalized.vcf.gz

# Filter for alleles that passed filtering and strip samples
echo `date`: Filtering VCF
cut -f 2-4 -d , $FILTERED_STR | tail -n+2 | sort | uniq | sed 's/:/\t/' | sed 's/,/\t/g' > $OUTPUT_DIR/Filtered_alleles.txt
$BCFTOOLS/bcftools view -R $OUTPUT_DIR/Filtered_alleles.txt -G $OUTPUT_DIR/Normalized.vcf.gz | $BCFTOOLS/bcftools sort | bgzip > $OUTPUT_DIR/Filtered.GangSTR.vcf.gz
tabix -p vcf $OUTPUT_DIR/Filtered.GangSTR.vcf.gz

# Use ANNOVAR to identify genic effects of STRs
echo `date`: Starting ANNOVAR
$ANNOVAR/table_annovar.pl \
	$OUTPUT_DIR/Filtered.GangSTR.vcf.gz \
	$ANNOVAR/humandb \
	-buildver hg38 \
	-out $OUTPUT_DIR/ANNOVAR \
	-remove -protocol wgEncodeGencodeBasicV38 \
	-operation g \
	-nastring . \
	-vcfinput -polish

bgzip $OUTPUT_DIR/ANNOVAR.hg38_multianno.vcf
tabix -p vcf $OUTPUT_DIR/ANNOVAR.hg38_multianno.vcf.gz

# Extract annotations per locus
echo `date`: Extracting annotations
$BCFTOOLS/bcftools query -HH -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/Func.wgEncodeGencodeBasicV38\t%INFO/Gene.wgEncodeGencodeBasicV38\t%INFO/ExonicFunc.wgEncodeGencodeBasicV38\n" $OUTPUT_DIR/ANNOVAR.hg38_multianno.vcf.gz > $OUTPUT_DIR/Parsed_ANNOVAR.txt

echo Ended at `date`
