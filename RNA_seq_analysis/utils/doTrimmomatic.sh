#trime adaptor 
#use slurm job

#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=do_trim_ipsc_npc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir /data6/johnathan/RNA_seq/ # TODO: set dir to cache dir
#SBATCH -o /data6/johnathan/RNA_seq/slurm/logs/out_trim_ipsc_npc_%a.log # TODO: set slurm output file
#SBATCH -e /data6/johnathan/RNA_seq/slurm/logs/err_trim_ipsc_npc%a.log # TODO: set slurm input file
#SBATCH --exclude=sarah
#SBATCH --exclude=laila
#SBATCH --array 1-6

export HOME="#path"

source /opt/anaconda/bin/activate #path, activate environment

echo $HOSTNAME

FILE=$(sed -n "$SLURM_ARRAY_TASK_ID"p #path_for_log_files)

MATE1="/path_for_fastq/${FILE}_R1_001.fastq.gz" #make a file for MATE1 and MATE2
MATE2="/path_for_fastq/${FILE}_R2_001.fastq.gz"
echo $MATE1
echo "======================"
echo $MATE2
echo "======================"
echo "======================"

out_R1_trimmed="/path_for_output/trimmomatic/trimmed/${FILE}_R1_trimmed.fastq.gz"
out_R1_unpaired="/path_for_output/unpaired/${FILE}_R1_unpaired.fastq.gz"
out_R2_trimmed="/path_for_output/trimmomatic/trimmed/${FILE}_R2_trimmed.fastq.gz"
out_R2_unpaired="/path_for_output/trimmomatic/unpaired/${FILE}_R2_unpaired.fastq.gz"
echo $out_R1_trimmed
echo "-------------------"
echo $out_R1_unpaired
echo "-------------------"
echo $out_R2_trimmed
echo "-------------------"
echo $out_R2_unpaired
echo "-------------------"
echo "-------------------"

trimmomatic PE -phred33 -threads 20 $MATE1 $MATE2 $out_R1_trimmed $out_R1_unpaired $out_R2_trimmed $out_R2_unpaired ILLUMINACLIP:/path_for_adaptor_sequenth/TruSeq3-PE.fa:2:30:10:7:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#adaptor sequence for this study
#>PrefixPE/1
#AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#>PrefixPE/2
#AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
