#slurmjob

#!/bin/bash
#SBATCH --account=girirajan # TODO: set account name
#SBATCH --partition=girirajan # TODO: set slurm partition
#SBATCH --job-name=redo_do_kallisto_ipsc_npc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=20G
#SBATCH --chdir  # TODO: set dir to cache dir
#SBATCH -o # TODO: set slurm output file
#SBATCH -e # TODO: set slurm input file
#SBATCH --exclude=sarah
#SBATCH --exclude=laila
#SBATCH --array 1-6

export HOME="/path_for_output"

source /opt/anaconda/bin/activate #path for kallisto environment 

echo $HOSTNAME

OUTDIR=$(sed -n "$SLURM_ARRAY_TASK_ID"p /slurm/files/ipsc_npc/2024-7-1_outdirs.txt)
MATE1=$(sed -n "$SLURM_ARRAY_TASK_ID"p /slurm/files/ipsc_npc/2024-7-1_mate1)
MATE2=$(sed -n "$SLURM_ARRAY_TASK_ID"p /slurm/files/ipsc_npc/2024-7-1_mate2) #mate 1 and mat2 are file names after trimmomatic

echo $MATE1
echo $MATE2

kallisto quant -i Homo_sapiens.GRCh38.cdna.all.kidx -o $OUTDIR -b 100 -t 5 $MATE1 $MATE2
