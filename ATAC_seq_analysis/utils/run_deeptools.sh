#FOXG1 site encode:ENCFF085GKB, ENCFF591AJU, ENCFF287YGM
#JUN site encode: ENCFF026KKS, ENCFF488BJY, ENCFF543DLZ, ENCFF550JRZ, ENCFF708RXW

#!/bin/bash
#SBATCH --job-name=npc_a38_cc
#SBATCH --output=final.out
#SBATCH --error=final.err
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --exclude ramona,durga,laila,qingyu

echo "Running job..."

#source /data5/jiawan/miniconda3/bin/activate
computeMatrix reference-point \
  --referencePoint center \
  -R /data5/jiawan/deeptool/foxg1_site.bed\
  -S /path_for_bigwig1.bw /path_for_bigwig2.bw /path_for_bigwig3.bw /path_for_bigwig4.bw \
  -o matrix.gz \
  --beforeRegionStartLength 2000 \
  --afterRegionStartLength 2000 \
  --sortRegions descend \
  --missingDataAsZero \
  --outFileNameMatrix matrix.txt

  echo "Job completed."

plotHeatmap -m matrix.gz -out heatmap.pdf  --colorMap Purples --refPointLabel "FOXG1 Motif"  --heatmapWidth 8 --heatmapHeight 16

plotProfile -m matrix.gz -out Profile.pdf --plotTitle "ATAC-seq Signal at FOXG1 Motifs" 
