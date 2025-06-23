#!/bin/bash
#SBATCH --account=johnathan # TODO: set account name
#SBATCH --partition=johnathan # TODO: set slurm partition
#SBATCH --job-name=filter_with_annovar
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=50G
#SBATCH --chdir /data7/ipsc/16p12_1_del/str_calls # TODO: set dir to project dir
#SBATCH -o /data7/ipsc/16p12_1_del/str_calls/slurm/logs/5_out_%a.log # TODO: set slurm output file
#SBATCH -e /data7/ipsc/16p12_1_del/str_calls/slurm/logs/5_err_%a.log # TODO: set slurm input file
#SBATCH --exclude=sarah # TODO: set nodelist
#SBATCH --array 1

export HOME="/data7/ipsc/16p12_1_del/str_calls"

while read LINE; do
    sample_id=`echo $LINE | cut -d" " -f1`
    vcf_files=`echo $LINE | cut -d " " -f2`
    echo "${sample_id}"
    outdir="/data7/ipsc/16p12_1_del/str_calls/data/filter_with_annovar/${sample_id}"
    mkdir -p $outdir

   table_annovar.pl ${vcf_files} humandb -buildver hg38 -out ${outdir}/${sample_id} -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput -polish
done < /data7/ipsc/16p12_1_del/str_calls/slurm/files/5_smap.txt
