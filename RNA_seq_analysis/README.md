Scipts of RNA-Seq analysis in the paper.
1. Trim adaptors use doTrimmomatic.sh
2. align reads and calculate counts use doKallisto.sh
3. tximport.R for input abundance files from Kallisto for the downsteam analysis.
4. DEseq2.R for analyzing differentially expressed genes (DEGs)
5. IsoformSwitchAnalyzeR.R for analyzing altered usage of isoforms.
6. Normalized the raw counts using rawcounts_to_TPMandFPKM.R
7. WGCNA.R for weighed gene-coexpression network analysis.
