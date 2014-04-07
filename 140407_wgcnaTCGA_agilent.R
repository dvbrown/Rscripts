getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
library(WGCNA)
options(stringsAsFactors=F)
list.files()

#################################### Data input and clustering visualisations ##############################################
data = read.delim('140110_agilentNoNulls.txt')

# We now transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(data))
# Check for genes with missing values and with no variance
gag = goodSamplesGenes(datExpr0, verbose=3)
gag$allOK

#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers. 
#We use the function flashClust that provides faster hierarchical clustering than the standard function hclust.
sampleTree = flashClust(dist(datExpr0), method = "average")
#pdf('./wgcna/140407_clusterSamples.pdf', paper='a4')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2, cex=0.6)
#dev.off()

# There are a small number of outliers but I will not trim them for the timebeing

##################################### Selecting the threshold ######################################################
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function. Use an unsigned network to begin with
#sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, blockSize=4000)
#save.image('140407_workspace.rdata')

# Plot the results of the thresholding function
#pdf('140407_thresholdingPlots.pdf', paper='a4')
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

##################################### Automatic network construction ######################################################

# You want to choose the power value based on where the curve flattens out
# This step is very slow
# net = blockwiseModules(datExpr0, power = 6,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "gbm_tcga",
#                        verbose = 3)
save.image('./wgcna/140407_networkBuilt.RData')