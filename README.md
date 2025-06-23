# iPSC_integrated_framework
Pipelines and scripts used to perform the bioinformatics analyses to build integrated framework for functional dissection of variable expressivity in genetic disease.


# Respository Organization

This repository contains scripts for perform the bioinformatics analyses in the paper. Files are organized into six main directories, each with additional READMEs. Brief description will be shown here. The  six main directories are:

# 1_Figure_code

Figure_code folder contains scripts to generate plots for the figures in the paper.

generate_heatmap.R for Figure1(c,d), Figure4(g), Figure6(e), Suppfig1(g), Suppfig4(c), Suppfig11(a).

generate_dotplot.R for Figure 1(e), Figure 2(b,c), Suppfig 1(e), Suppfig 3(b), Suppfig 8(d) and Suppfig 11 (a)

barplot.R for Figure 3, Figure 4, Figure 6, Figure 5, Suppfig 8 and Suppfig 10.

volcanoplot.R for Suppfig 1c, Suppfig 2a and Suppfig 7d.

upsetplot.R for Figure 4 (b).

waterfallplot.R for Figure 4 (g), Suppfig 1 (f) and Suppfig 6(c).

densityplot. R for Figure 5 (b) and Suppfig 9 (a).

lineplot.py for Suppfig 2(c).

violinplot.R for Figure 3 (a), Suppfig 4(a), Suppfig 5

# 2_Enrichment_analysis

Scripts for enrichment analysis in the paper.

Enrichment_disgenet.R for Suppfig 1 (e) and Suppfig 3(b).

Enrichment_in_published_datases.R for Suppfig 1 (e) and Suppfig 3(b).

GSEA.R for gene-set enrichment analysis in the paper. 

ORA.R for over-representation analysis in the paper.

# 3_RNA_seq_analysis

Trim adaptors use doTrimmomatic.sh.

align reads and calculate counts use doKallisto.sh.

tximport.R for input abundance files from Kallisto for the downsteam analysis.

DEseq2.R for analyzing differentially expressed genes (DEGs).

IsoformSwitchAnalyzeR.R for analyzing altered usage of isoforms.

Normalized the raw counts using rawcounts_to_TPMandFPKM.R.

WGCNA.R for weighed gene-coexpression network analysis.

# 4_ATAC_seq_analysis

Encode ATAC-seq pipline (https://github.com/ENCODE-DCC/atac-seq-pipeline) for preproceesing.

makeATACTable.pl for summarize counts per regions.

DEseq2_ATAC.R for analyzing diffPeaks

run_mea_homer.sh for motif enrichment analysis

get_nearest_gene.py, chipseeker (use galaxy psu) and ABC enhancer-gene prediction (https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) for annotate peaks.

run_deeptools.sh for normalize binding site activity.

# 5_WGS_analysis

call_CNV for calling, filtering and annotating CNVs.

call_SNV for calling, filtering and annotation SNVs.

call_STR for calling, filtering and annotation STRs.

Fisher_exact_test for overlap analysis between WGS variants and DEGs or Diffpeaks.

Gene_interaction_analysis for finding variants and genes with additional changes in gene expression and chromatin accessibility in families.

# 6_linear_regression

Scripts for linear regression analysis in Figure 7.

Input tables and data are in Supplementary table 10.

Modified from https://github.com/girirajanlab/16p12.1-Deletion-WGS
