getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
library(WGCNA)
library(RColorBrewer)
options(stringsAsFactors=F)
list.files()

#################################### Data input and clustering visualisations ##############################################
data = read.delim('140110_agilentNoNulls.txt')

# We now transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(data))
dat = as.matrix(datExpr0)
save.image('140415_justTheAgilentData.RData')
########################################################## Standard correlations ###############################################
load('140415_justTheAgilentData.RData')

# Calculate the correlation between PROM1 expression and all the genes in TCGA GBM
prom1CorrPval = corAndPvalue(x=datExpr0[,'PROM1'], y=dat)

#Extract the correlation and p-value from the returned list
prom1C = prom1CorrPval$cor
prom1Cpower = prom1C^2

prom1P = prom1CorrPval$p
prom1FDR = p.adjust(prom1P, method='fdr')

# par(mfrow=c(2,1))
# hist(prom1Cpower, main='Prom1 correlations', breaks='FD', xlab='Weighted correlation values')
# hist(prom1FDR, main='Prom1 p-values', breaks='FD', xlab='FDR corrected p-values')
# par(mfrow=c(1,1))

result = t(rbind(prom1C, prom1Cpower, prom1P, prom1FDR))
colnames(result) = c('correlation', 'weighted_correlation', 'p-value', 'FDR')
#write.table(result, './wgcna/140415_standardCorrelation.txt', sep='\t')

par(mfrow=c(2,2))
hist(prom1C, main='Prom1 correlations', breaks='FD', xlab='Raw correlation values')
hist(prom1Cpower, main='Prom1 weighted correlations', breaks='FD', xlab='R squared values')
hist(prom1P, main='Prom1 p-values', breaks='FD', xlab='Raw p-values')
hist(prom1FDR, main='Prom1FDR corrected p-values', breaks='FD', xlab='FDR corrected p-values')
par(mfrow=c(1,1))

qqnorm(result[,1], main='Distrubution of PROM1 correlated genes')
qqline(result[,1])

######################################### Build a heatmap of the correlated genes using the square matrix #####################

# Subset the dataframe with correlation values for those with high correlation and significance
prom1Cgenes = result[result[,2] > 0.1 & result[,4] < 0.05,]
length(row.names(prom1Cgenes))

prom1CgenesNames = row.names(prom1Cgenes)
prom1CgenesNames

# Build the network adjacency
# Use the top correlated genes with PROM1 and measure their correlation with the transcriptome
adjacencyProm1 = adjacency(datExpr0, 
                            selectCols = prom1CgenesNames, #for correlation networks only (see below); can be used to select genes whose adjacencies will be calculated. Should be either a numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are to be used.
                            type = "unsigned", power = 6, corFnc = "cor", #corOptions = "use = 'p'",
                            distFnc = "dist", distOptions = "method = 'euclidean'")

# Make the adjacency matrix square
squareAdjacency = (adjacencyProm1[colnames(adjacencyProm1),])

#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix,
#and calculate the corresponding dissimilarity:

# Calculate the topological overlap matrix
similarity = TOMsimilarity(squareAdjacency, TOMType='unsigned', verbose=3)
row.names(similarity) = row.names(squareAdjacency)
colnames(similarity) = row.names(squareAdjacency)

cc = brewer.pal(9, 'YlOrRd')
#heatmap(log10(squareAdjacency), main='Most adjacent genes to CD133 (n = 134)', Rowv=NA, sym=TRUE)
heatmap(log10(similarity), main='Most similar genes to CD133 (n = 134)', Rowv=NA, sym=TRUE, col=cc, cexRow=0.5, cexCol=0.5)
###########################################################################################################################

######################################### Visulaising the network that was built ##########################################


dissTOM = 1-similarity
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#Branches of the dendrogram group together densely interconnected, highly co-expressed genes. Module identifcation amounts to the identification of individual branches
#There are several methods for branch cutting; our standard method is the Dynamic Tree Cut from the package dynamicTreeCut

# Module identification using dynamic tree cut. This is the most basic method and returns 3 modules when the cutHeight is 0.999 (default 0.99)
dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=0.999)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")