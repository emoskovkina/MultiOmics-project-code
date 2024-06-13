##########################################################################################
##########################################################################################
######################### APPLIED COMPUTATIONAL MULTI-OMICS ##############################
################################ TUTORIAL 4 - RNAseq #####################################
##########################################################################################
##########################################################################################


#This practical tutorial will perform differentia gene expression analysis using R packages. 
#Most steps will be performed using edgeR R package, which complete tutorial can be found here: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


#########################################################################################################################
# Install packages from R and BioConductor repositories ####
#########################################################################################################################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("tximport","rhdf5", "PCAtools", "edgeR"))

install.packages(c("gplots", "RColorBrewer", "openxlsx"))


#########################################################################################################################
# Import read counts to R ####
# We will import the read counts obtained with Kallisto (Turorial 3) using the tximport R package
#########################################################################################################################

library(tximport)

# First, set the R working directory to the folder containing Kallisto results 
setwd("~/nova/mo") #Complete with the path in your computer

# Load the file containing the EnsemblIDs and gene names
load('gencode.v38_geneInfo.RData');
head(geneInfo)

# Define filenames to import
setwd("~/nova/mo/project/kallisto_output")
folders <- dir(pattern="_rep")
fileNames <- file.path(".", folders, 'abundance.tsv');
names(fileNames) <- folders 

# Import all data using tximport function
readCounts_gene <- tximport(fileNames, type = 'kallisto',txOut = F,tx2gene=geneInfo[,1:2], ignoreTxVersion = T,countsFromAbundance = "lengthScaledTPM")

# Explore the object created with tximport
names(readCounts_gene)
View(head(readCounts_gene$abundance)) # matrix containing TPMs
View(head(readCounts_gene$counts)) # matrix containing read counts
head(readCounts_gene$length) 
head(readCounts_gene$countsFromAbundance) # indicates the type of normalization applied to get the values in "abundance"

# To facilitate interpretation of downstream results, replace EnsemblIDs with Official gene symbol (i.s gene name)
geneNames <- geneInfo[match(rownames(readCounts_gene$abundance), geneInfo$Ensembl_GeneID),"GeneSymbol"]
rownames(readCounts_gene$abundance) <- rownames(readCounts_gene$counts) <- geneNames
head(readCounts_gene$abundance)

head(readCounts_gene$counts)


# Save readCounts in an RData file for downstream analyses
save(readCounts_gene, file="readCounts_gene.RData")

#########################################################################################################################
# Explore transcriptome profiles using PCA ####
#########################################################################################################################

library(PCAtools)

# Prepare gene expression matrix for PCA analysis: Step 1) get log TPMs (we add 1 TPM to each gene to avoid infinite values after log)
logTPMs <- log2(readCounts_gene$abundance+1)
head(readCounts_gene$abundance)
boxplot(logTPMs)
# Prepare gene expression matrix for PCA analysis: Step 2) remove duplicated genes
uniqueGenes <- unique(rownames(logTPMs))
logTPMs <- logTPMs[uniqueGenes,]
boxplot(logTPMs)

#Prepare metadata with sample type
colnames(logTPMs)
sampleTypes <- gsub("_rep[123]", "", colnames(logTPMs))
sampleTypes
metaData <- data.frame(sampleTypes); rownames(metaData) <- colnames(logTPMs)
metaData


#Run PCA
pca.res <- pca(logTPMs, metadata=metaData)

#Plot variance explained by each component
screeplot(pca.res)

#Plot 2 selected components/eigenvectors
biplot(pca.res)
biplot(pca.res, colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res, lab="", colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res, x="PC1", y="PC3",lab="",colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3


#Plot several components
pairsplot(pca.res, colby="sampleTypes")
pairsplot(pca.res)

# Plot the component loadings and label genes most responsible for variation
# NOTE: Loadings are interpreted as the coefficients of the linear combination of the initial variables from which the principal components are constructed.
# From a statistical point of view, the loadings are equal to the coordinates of the variables divided by the square root of the eigenvalue associated with the component.
PC_genes <- pca.res$loadings
PC1_genes <- PC_genes[order(PC_genes$PC1, decreasing=T),]
head(PC1_genes)
tail(PC1_genes)
plotloadings(pca.res, components = c("PC1", "PC2", "PC3"),rangeRetain =0.1) #retaining 1% of the loadings per PC

# Produce barplot to confirm expression levels of genes associated with Principal Component 1
plotCol <- rep(c("steelblue", "gray", "orange", "darkviolet"), each=3)
colnames(logTPMs)
barplot(logTPMs["WWTR1-AS1",], col=plotCol, las=2, main="WWTR1-AS1", ylab="Expression levels (logTPMs)") #Gene positively correlated with PC1
legend("topleft", fill=unique(plotCol), legend=c("Control", "TAZ_S89A_overexpression","YAP_S127A_overexpression","YAP_S94A_S127A_overexpression"), bty="n")
barplot(logTPMs["FABP4",], col=plotCol, las=2, main="FABP4", ylab="Expression levels (logTPMs)") #Gene negatively correlated with PC1
legend("topright", fill=unique(plotCol), legend=c("Control", "TAZ_S89A_overexpression","YAP_S127A_overexpression","YAP_S94A_S127A_overexpression"), bty="n")


# You can save some plot in a pdf file
pdf("PCA_plots.pdf") # open pdf file
biplot(pca.res, colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
plotloadings(pca.res, components = c("PC1", "PC2", "PC3"),rangeRetain =0.1) 
dev.off() # closes pdf file


#########################################################################################################################
# Differential expression analysis using edge R package ####
#########################################################################################################################


library(edgeR)

# Create DGEList data class (specific for edgeR package)
y <- DGEList(counts=readCounts_gene$counts, group= sampleTypes)
row.names(readCounts_gene$counts)

# Filter out lowly expressed genes
keep <- filterByExpr(y, group=sampleTypes)
y <- y[keep, ,keep.lib.sizes=FALSE]

# Normalize for library sizes
y <- calcNormFactors(y)
y$samples


# Define the design and contrast matrix based on the experimental design, meaning define which comparison to be made
# See more detailed information here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
design_matrix <- model.matrix(~0+sampleTypes)
colnames(design_matrix) <- gsub("sampleTypes", "", colnames(design_matrix))
rownames(design_matrix) <- colnames(logTPMs)
design_matrix

############################

contrast_matrix <- makeContrasts(TAZ_S89A_overexpression-Control, levels = design_matrix) # NOTE: REPPLACE HERE WITH THE NAMES OF YOUR SAMPLES
contrast_matrix


# Estimate the dispersion (Biological coefficient of variation) and Fit model
# edgeR uses the negative binomial (NB) distribution to model the read counts for each gene in each sample. The dispersion parameter of the NB distribution accounts for variability between biological replicates (McCarthy, Chen, and Smyth 2012). edgeR estimates an empirical Bayes moderated dispersion for each individual gene. It also estimates a common dispersion, which is a global dispersion estimate averaged over all genes, and a trended dispersion where the dispersion of a gene is predicted from its abundance.
# The vertical axis of the plotBCV plot shows square-root dispersion, also known as biological coefficient of variation (BCV) (McCarthy, Chen, and Smyth 2012).
# For RNA-seq studies, the NB dispersions tend to be higher for genes with very low counts. The dispersion trend tends to decrease smoothly with abundance and to asymptotic to a constant value for genes with larger counts. From our past experience, the asymptotic value for the BCV tends to be in range from 0.05 to 0.2 for genetically identical mice or cell lines, whereas somewhat larger values (>0.3) are observed for human subjects.
y <- estimateDisp(y, design_matrix)
plotBCV(y)
fit <- glmQLFit(y,design_matrix) # Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 


# Differential expression analysis
qlf <- glmQLFTest(fit, contrast =  contrast_matrix) # test for differential expression between the experimental groups using quasi-likelihood F-test

# Get differentially expressed genes (DEGs)
summary(decideTests(qlf,p.value = 0.05,adjust.method = "fdr"))
DEGs <- topTags(qlf, nrow(qlf$table), p.value=0.05, adjust.method = "fdr")$table
upReg <- rownames(DEGs)[which(DEGs$logFC > 0)]
downReg <- rownames(DEGs)[which(DEGs$logFC < 0)]

# Write table with DEGs
library(openxlsx)
write.xlsx(DEGs, file = "DEGs_TAZ_S89A_overexpression.xlsx", rowNames=T)
row.names(DEGs)
write.table(row.names(DEGs), file="geneNames_TAZ_S89A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(upReg, file="geneNames_upReg_TAZ_S89A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(downReg, file="geneNames_downReg_TAZ_S89A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(data.frame(row.names(DEGs),DEGs$logFC), file="rnk_TAZ_S89A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")

DEGs_with_cutoff <- rownames(DEGs)[which(DEGs$logFC > 1 | DEGs$logFC < -1)]
length(DEGs_with_cutoff)
write.table(DEGs_with_cutoff, file="geneNames_TAZ_S89A_overexpression_cutoff_logFC_1.txt", quote=F, row.names=F, col.names=F, sep="\t")
upReg_c <- rownames(DEGs)[which(DEGs$logFC > 1)]
downReg_c <- rownames(DEGs)[which(DEGs$logFC < 1)]
write.table(upReg_c, file="geneNames_upReg_TAZ_S89A_overexpression_cutoff_logFC_1.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(downReg_c, file="geneNames_downReg_TAZ_S89A_overexpression_cutoff_logFC_1.txt", quote=F, row.names=F, col.names=F, sep="\t")

# Produce Volcano Plot (customized)
allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table

head(allGenes)
write.table(row.names(allGenes), file="geneNames.txt", quote=F, row.names=F, col.names=F, sep="\t")
plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
abline(h=-log10(0.05), v=c(-1,1),lty=2)
points(plotData[upReg,], col="tomato", pch=20)
points(plotData[downReg,], col="steelblue", pch=20)
text(plotData[upReg[1:5],], labels=upReg[1:5],col="tomato", pos=sample(c(1:3), size=10, replace=T), cex=0.8)
text(plotData[downReg[1:5],], labels=downReg[1:5],col="steelblue", pos=sample(c(1:2), size=10, replace=T), cex=0.8)

# Produce MA plot (customized)
plotData <- cbind(allGenes$logCPM, allGenes$logFC); rownames(plotData) <- rownames(allGenes)
plot(plotData, pch=20,col="gray",xlab="Mean Expression Levels (log CPMs)", ylab="Biological Variation (log2 Fold-Change)")
abline(h=c(-1,1),lty=2)
points(plotData[upReg,], col="tomato", pch=20)
points(plotData[downReg,], col="steelblue", pch=20)
text(plotData[upReg[1:5],], labels=upReg[1:5],col="tomato", pos=sample(c(1:3), size=10, replace=T), cex=0.8)
text(plotData[downReg[1:5],], labels=downReg[1:5],col="steelblue", pos=sample(c(1:2), size=10, replace=T), cex=0.8)


# Heatmap with Top DEGs
library(gplots); library(RColorBrewer)
plotCol_exp <- brewer.pal(9, "Purples")
heatmap.2(logTPMs[c(upReg[1:10], downReg[1:10]),], scale="row", trace="none", density.info="none", ColSideColors = plotCol, col=plotCol_exp)


###################################################3

contrast_matrix <- makeContrasts(YAP_S94A_S127A_overexpression-Control, levels = design_matrix) # NOTE: REPPLACE HERE WITH THE NAMES OF YOUR SAMPLES
contrast_matrix


# Differential expression analysis
qlf <- glmQLFTest(fit, contrast =  contrast_matrix) # test for differential expression between the experimental groups using quasi-likelihood F-test

# Get differentially expressed genes (DEGs)
summary(decideTests(qlf,p.value = 0.05,adjust.method = "fdr"))
DEGs <- topTags(qlf, nrow(qlf$table), p.value=0.05, adjust.method = "fdr")$table
upReg <- rownames(DEGs)[which(DEGs$logFC > 0)]
downReg <- rownames(DEGs)[which(DEGs$logFC < 0)]

# Write table with DEGs
write.xlsx(DEGs, file = "DEGs_YAP_S94A_S127A_overexpression.xlsx", rowNames=T)
row.names(DEGs)
write.table(row.names(DEGs), file="geneNames_YAP_S94A_S127A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(upReg, file="geneNames_upReg_YAP_S94A_S127A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(downReg, file="geneNames_downReg_YAP_S94A_S127A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")
write.table(data.frame(row.names(DEGs),DEGs$logFC), file="rnk_YAP_S94A_S127A_overexpression.txt", quote=F, row.names=F, col.names=F, sep="\t")

# Produce Volcano Plot (customized)
allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table

head(allGenes)
write.table(row.names(allGenes), file="allGeneNames_YAP_S94A_S127A.txt", quote=F, row.names=F, col.names=F, sep="\t")
plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
abline(h=-log10(0.05), v=c(-1,1),lty=2)
points(plotData[upReg,], col="tomato", pch=20)
points(plotData[downReg,], col="steelblue", pch=20)
text(plotData[upReg[1:5],], labels=upReg[1:5],col="tomato", pos=sample(c(1:3), size=10, replace=T), cex=0.8)
text(plotData[downReg[1:5],], labels=downReg[1:5],col="steelblue", pos=sample(c(1:2), size=10, replace=T), cex=0.8)

# Produce MA plot (customized)
plotData <- cbind(allGenes$logCPM, allGenes$logFC); rownames(plotData) <- rownames(allGenes)
plot(plotData, pch=20,col="gray",xlab="Mean Expression Levels (log CPMs)", ylab="Biological Variation (log2 Fold-Change)")
abline(h=c(-1,1),lty=2)
points(plotData[upReg,], col="tomato", pch=20)
points(plotData[downReg,], col="steelblue", pch=20)
text(plotData[upReg[1:5],], labels=upReg[1:5],col="tomato", pos=sample(c(1:3), size=10, replace=T), cex=0.8)
text(plotData[downReg[1:5],], labels=downReg[1:5],col="steelblue", pos=sample(c(1:2), size=10, replace=T), cex=0.8)


# Heatmap with Top DEGs
library(gplots); library(RColorBrewer)
plotCol_exp <- brewer.pal(9, "Purples")
heatmap.2(logTPMs[c(upReg[1:10], downReg[1:10]),], scale="row", trace="none", density.info="none", ColSideColors = plotCol, col=plotCol_exp)


