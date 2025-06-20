library(WGCNA)
library(FactoMineR) 
library(factoextra) 
library(tidyverse) 
library(data.table)
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores())

tpm <- family_npc_counts_tpm #input tpm
datTraits <- datTraits_npc  # input sample information 

data <- log2(tpm+1)
keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:10000],]

datExpr0 <- as.data.frame(t(keep_data))

#quality_check
sampleTree <- hclust(dist(datExpr0), method = "average")
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
cex.axis = 1, cex.main = 1,cex.lab=1)
sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)),
colors = rainbow(length(table(datTraits$group))),
signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
pdf("Sample dendrogram and trait.pdf",width = 8,height = 6)
plotDendroAndColors(sampleTree, sample_colors,
groupLabels = "trait",
cex.dendroLabels = 0.8,
marAll = c(1, 4, 3, 1),
cex.rowText = 0.01,
main = "Sample dendrogram and trait" )
dev.off()

group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F)
pca <- fviz_pca_ind(dat.pca,
title = "Principal Component Analysis",
legend.title = "Groups",
geom.ind = c("point","text"), #"point","text"
pointsize = 2,
labelsize = 4,
repel = TRUE, 
col.ind = group_list, 
axes.linetype=NA, # remove axeslines
mean.point=F
) +
theme(legend.position = "none")+ # "none" REMOVE legend
coord_fixed(ratio = 1)
pca
ggsave(pca,filename= "Sample PCA analysis.pdf", width = 8, height = 8)


datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#pick_power_for_unsigend

R.sq_cutoff = 0.8 
if(T){
# Call the network topology analysis function

powers <- c(seq(1,20,by = 1), seq(22,30,by = 2))
sft <- pickSoftThreshold(datExpr,
networkType = "unsigned",# can also be signed or signed hybrid
powerVector = powers,
RsquaredCut = R.sq_cutoff,
verbose = 5)
#SFT.R.sq > 0.8 , slope ≈ -1
pdf("power-value_unsigned.pdf",width = 16,height = 12)

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=R.sq_cutoff ,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=100,col="red")
dev.off()
}

power = sft$powerEstimate


#build_network
allowWGCNAThreads()
net <- blockwiseModules(
datExpr,
power = power,
maxBlockSize = ncol(datExpr),
corType = "pearson", 
networkType = "unsigned", #use signed or signed hybrid with the setting when picking power.
TOMType = "unsigned",
minModuleSize = 100, 
mergeCutHeight = 0.25, 
numericLabels = TRUE,
saveTOMs = F,
verbose = 3
)
table(net$colors)

# Convert labels to colors for plotting
moduleColors <- labels2colors(net$colors)
table(moduleColors)
# Plot the dendrogram and the module colors underneath
pdf("genes-modules_ClusterDendrogram_unsigned.pdf",width = 16,height = 12)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

gene_module <- data.frame(gene=colnames(datExpr),
module=moduleColors)
write.csv(gene_module,file = "gene_moduleColors_unsigned.csv",row.names = F)
MES0 <- moduleEigengenes(datExpr, moduleColors)$eigengnes #calculate the eigengene score
MEs_unsigned <- orderMEs(MES0) #put close eigenvectors next to each other
write.csv(MEs_unsigned,file = "MEs_unsigned.csv")

#draw heatmap using MEs such as Suppfig 6b
datTraits$group <- as.factor(datTraits$group)
design <- model.matrix(~0+datTraits$group)
colnames(design) <- levels(datTraits$group) #get the group
MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes #Calculate module eigengenes.
MEs <- orderMEs(MES0) #Put close eigenvectors next to each other
moduleTraitCor <- cor(MEs,design,use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
signif(moduleTraitPvalue,1),")")
dim(textMatrix) <- dim(moduleTraitCor)
pdf("step4_Module-trait-relationship_heatmap.pdf",
width = 2*length(colnames(design)),
height = 0.6*length(names(MEs)) )
par(mar=c(5, 9, 3, 3)) 
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = colnames(design),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = F,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = F,
cex.text = 0.5,
zlim = c(-1,1),
main = "Module-trait relationships")
dev.off()

#output files for cytoscape such as Suppfig 11b
gene <- colnames(datExpr)
inModule <- moduleColors==module #input module color of interest
modgene <- gene[inModule]

TOM <- TOMsimilarityFromExpr(datExpr,power=power)
modTOM <- TOM[inModule,inModule]
dimnames(modTOM) <- list(modgene,modgene)

nTop = 100
IMConn = softConnectivity(datExpr[, modgene]) #calculate connectivity 
top = (rank(-IMConn) <= nTop) #filter for the top
filter_modTOM <- modTOM[top, top]

cyt <- exportNetworkToCytoscape(filter_modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.15, #can be changed 
nodeNames = modgene[top],
nodeAttr = moduleColors[inModule][top])

