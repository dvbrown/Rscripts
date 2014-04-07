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

##################################### Automatic network construction ######################################################
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)