getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
library(WGCNA)
options(stringsAsFactors=F)
list.files()

#################################### Data input and clustering visualisations ##############################################
data = read.delim('140110_agilentNoNulls.txt')

# We now transpose the expression data for further analysis.
datExpr = as.data.frame(t(data))

# Use only 100 patients to make things run quicker. CHANGE THIS to full set when optimised
datExpr0 = datExpr[c(1:100),]
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

# Choose a height cut that will remove the oending sample, say 140 (the red line in the plot), and use a branch cut at that height.
# Plot a line to show the cut
abline(h = 140, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

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
#save.image('./wgcna/140407_networkBuilt.RData')
load('./wgcna/140407_networkBuilt.RData')

# Identify how many modules there are and how bif theu are.
table(net$colors)

# View the dendogram you can cut the tree without recutting it.
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# We now save the module assignment and module eigengene information necessary for subsequent analysis.
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

# That wasn't particularly informative and there are many modules

##################################### Manual network construction ######################################################

# We now calculate the adjacencies, using the soft thresholding power 6:
softPower = 6
adjacency = adjacency(datExpr0, power = softPower)