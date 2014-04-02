getwd()
setwd('~/Documents/public-datasets/WGCNA_tutorial/FemaleLiver-Data')
library(WGCNA)
options(stringsAsFactors=F)

# Let R do some multithreading
allowWGCNAThreads()
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv");
# Take a quick look at what is in the data set:
dim(femData)
names(femData)

# We now remove the auxiliary data and transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH;
# Remove the metaData columns
rownames(datExpr0) = names(femData)[-c(1:8)];

# Check for genes with missing values and with no variance
gag = goodSamplesGenes(datExpr0, verbose=3)
gag$allOK

if (!gsg$allOK)
{
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers. We use the function flashClust that provides faster hierarchical clustering than the standard function hclust.
sampleTree = flashClust(dist(datExpr0), method = "average")
# Plot the sample tree
pdf('./140402_clucterSamples.pdf', paper='a4')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# It appears there is one outlier (sample F2_221, see Fig. 1). One can remove it by hand, or use an automatic approach.
# Choose a height cut that will remove the oending sample, say 15 (the red line in the plot), and use a branch cut at that height.
# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# The variable datExpr now contains the expression data ready for network analysis.
# Read in clinical data
traitData = read.csv("ClinicalTraits.csv");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Mice);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbagee();

# We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable
# datTraits. Before we continue with network construction and module detection, we visualize how the clinical traits
# relate to the sample dendrogram.