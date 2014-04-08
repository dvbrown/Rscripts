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

You want to choose the power value based on where the curve flattens out
This step is very slow
net = blockwiseModules(datExpr0, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "gbm_tcga",
                       verbose = 3)
# MEs = a data frame containing module eigengenes of the found modules (given by colors).

save.image('./wgcna/140407_networkBuilt.RData')
load('./wgcna/140407_networkBuilt.RData')

# Identify how many modules there are and how big theu are.
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

##################################### Identify which genes are highly correlated with modules ######################################################
dat = as.matrix(datExpr0)
geneModuleMembership = as.data.frame(corFast(dat, MEs, use = "p"))
geneModuleMembership = abs(geneModuleMembership)

# Subset the module membership for CD133
prom1 = t(geneModuleMembership['PROM1',])
# It appears to belong to at least 2 modules

prom1 = sort(prom1[,1], decreasing=TRUE)

prom1Members = geneModuleMembership$ME19
names(prom1Members) = row.names(geneModuleMembership)
hist(prom1Members)

##################################### Try calculate network adjaceny again ######################################################
# Find the index of prom1 using tail. Need to find a computational way
prom1 = 12404
cnames = colnames(datExpr0)
# This works as a method to subset the gene expression matrix
prom1 = ifelse(cnames %in% 'PROM1', TRUE, FALSE)
egfr = ifelse(cnames %in% 'EGFR', TRUE, FALSE)

# Build a vector that matches the presence of an index

# Calculate the adjacency for only prom1.  (correlation or distance) network adjacency
adj = adjacency(datExpr0, 
                selectCols = prom1, #for correlation networks only (see below); can be used to select genes whose adjacencies will be calculated. Should be either a numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are to be used.
                type = "unsigned", power = 6, corFnc = "cor", #corOptions = "use = 'p'",
                distFnc = "dist", distOptions = "method = 'euclidean'")

hist(adj, breaks='FD', ylim=c(0,20), main='Frequency of adjacency values', xlab='Adjacency')

prom1Genes = adj[adj > 0.5]

# adjE = adjacency(datExpr0, 
#                 selectCols = egfr, #for correlation networks only (see below); can be used to select genes whose adjacencies will be calculated. Should be either a numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are to be used.
#                 type = "unsigned", power = 6, corFnc = "cor", #corOptions = "use = 'p'",
#                 distFnc = "dist", distOptions = "method = 'euclidean'")


##################################### Start here tomorrow ######################################################
adj = adjacency(datExpr0, 
                selectCols = null, #for correlation networks only (see below); can be used to select genes whose adjacencies will be calculated. Should be either a numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are to be used.
                type = "unsigned", power = 6, corFnc = "cor", #corOptions = "use = 'p'",
                distFnc = "dist", distOptions = "method = 'euclidean'")

TOM = TOMsimilarity(adj);
dissTOM = 1-TOM
