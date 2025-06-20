#slurm job
#!/bin/bash
#SBATCH --job-name=ipsc_3110_331_321_342_up
#SBATCH --output=my_job_ipsc_3110_331_321_342_up.out
#SBATCH --error=my_job_ipsc_3110_331_321_342_up.err
#SBATCH --time=300:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

echo "Running job..."

#source /opt/anaconda/bin/activate /data5/deepro/miniconda3/envs/starrseq_mea/

findMotifsGenome.pl /path_for_dp_results_up.bed /GRCh38.p14.genome.fa /path_for_output -p 32 -size given -bg /path_for_atac_peaks_background.bed
echo "Job completed."
