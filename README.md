# iPSC_integrated_framework
Pipelines and scripts used to perform the bioinformatics analyses to build integrated framework for functional dissection of variable expressivity in genetic disease.


# Respository Organization

This repository contains scripts for perform the bioinformatics analyses in the paper. Files are organized into six main directories, each with additional READMEs. Brief description will be shown here. The  six main directories are:

# 1_Figure_code

Figure_code folder contains scripts to generate plots for the figures in the paper.

generate_heatmap.R for Figure1(c), Figure4(f-g), Figure6(e), Suppfig3(c), Suppfig5(d).

generate_dotplot.R for Figure 2(b,d), Figure 3(g), Suppfig 3(a-b,d), Suppfig 4(a-e), Suppfig 7(b) and Suppfig 9(d)

barplot.R for Figure 3, Figure 4, Figure 6, Suppfig 6, Suppfig 8, Suppfig 9 and Suppfig 11.

volcanoplot.R for Suppfig 1(d), Suppfig 2(a), Suppfig 5(c) and Suppfig 8(b).

upsetplot.R for Figure 4(b).

densityplot. R for Figure 5(a) and Suppfig 10(a).

lineplot.py for Suppfig 2(c).

violinplot.R for Figure 3, Suppfig 5, Suppfig 6.

# 2_Enrichment_analysis

Scripts for enrichment analysis in the paper.

Enrichment_disgenet.R for Suppfig 3(b) and Suppfig 4(b,d).

Enrichment_in_published_datases.R for Suppfig 3(a) and Suppfig 4(a,c).

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

Fisher_exact_test.py for Example code of Fisher's exact tests.

Gene_interaction_analysis for finding variants and genes with additional changes in gene expression and chromatin accessibility in families.


