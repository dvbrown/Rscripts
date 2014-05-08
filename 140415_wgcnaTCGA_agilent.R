getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/wgcna/')
library(WGCNA)
options(stringsAsFactors=F)
list.files()

#################################### Data input and clustering visualisations ##############################################

# data = read.delim('../140110_agilentNoNulls.txt')
# 
# # We now transpose the expression data for further analysis.
# datExpr0 = as.data.frame(t(data))
# dat = as.matrix(datExpr0)
# save.image('140415_justTheAgilentData.RData')

########################################################## Standard correlations ###############################################

load('140415_justTheAgilentData.RData')

# Calculate the correlation between PROM1 expression and all the genes in TCGA GBM
prom1CorrPval = corAndPvalue(x=dat[,'PROM1'], y=dat)

#Extract the correlation and p-value from the returned list
prom1C = prom1CorrPval$cor
prom1Cpower = prom1C^2

prom1P = prom1CorrPval$p
prom1FDR = p.adjust(prom1P, method='fdr')

par(mfrow=c(2,2))
hist(prom1C, main='Prom1 correlations', breaks='FD', xlab='Weighted correlation values')
hist(prom1Cpower, main='Prom1 person regression', breaks='FD', xlab='Weighted correlation values')
hist(prom1FDR, main='Prom1 p-values', breaks='FD', xlab='FDR corrected p-values')
par(mfrow=c(1,1))

result = t(rbind(prom1C, prom1Cpower, prom1P, prom1FDR))
colnames(result) = c('correlation', 'weighted_correlation', 'p-value', 'FDR')
#write.table(result, './wgcna/140415_standardCorrelation.txt', sep='\t')

# par(mfrow=c(2,2))
# hist(prom1C, main='Prom1 correlations', breaks='FD', xlab='Raw correlation values')
# hist(prom1Cpower, main='Prom1 weighted correlations', breaks='FD', xlab='R squared values')
# hist(prom1P, main='Prom1 p-values', breaks='FD', xlab='Raw p-values')
# hist(prom1FDR, main='Prom1FDR corrected p-values', breaks='FD', xlab='FDR corrected p-values')
# par(mfrow=c(1,1))
# 
# qqnorm(result[,1], main='Distrubution of PROM1 correlated genes')
# qqline(result[,1])

######################################### Subsample the data matrix to validate #####################


# Look up http://www.pmc.ucsc.edu/~mclapham/Rtips/resampling.htm

head(sample(dat))
# Subsample the dataframe taking only 100 patients
subDat = dat[sample(1:nrow(dat),100),] #randomizes a data frame; add replace=T for bootstrapping and n for subsampling
subsamplingCorr = corAndPvalue(x=subDat[,'PROM1'], y=subDat)

par(mfrow=c(1,2))
hist(prom1C, main='Prom1 correlations', breaks='FD', xlab='Weighted correlation values')
hist(subsamplingCorr$cor, main='Subsampled Prom1 correlations ', breaks='FD', xlab='Weighted correlation values')

hist(subsamplingCorr$cor^2, main='Subsampled Prom1 \nPearson regression', breaks='FD', xlab='Weighted correlation values')
hist(p.adjust(subsamplingCorr$p, 'fdr'), main='Subsampled Prom1 corrected p-values', breaks='FD', xlab='FDR corrected p-values')
par(mfrow=c(1,1))

subSampledResult = t(rbind(subsamplingCorr$cor, subsamplingCorr$cor^2, subsamplingCorr$p, p.adjust(subsamplingCorr$p, 'fdr')))

# Subset the dataframe with correlation values for those with high correlation and significance
subSampledProm1Genes = subSampledResult[subSampledResult[,2] > 0.1 & subSampledResult[,4] < 0.05,]

######################################### Visulaise the most correlated genes #####################
# par(mfrow=c(2,1))
# plot.default(dat[,'PROM1'], dat[,'PHF15'], ylab='PHF15', xlab='CD133', main='Genes coexpressed with CD133')
# correlations = lm(dat[,'PHF15'] ~ dat[,'PROM1'])
# abline(correlations, col='red')
# summary(correlations)
# text(locator(1), 'R-squared: 0.2175')
# 
# # Visualise the most next correlated gene
# plot.default(dat[,'PROM1'], dat[,'SOX11'], ylab='SOX11', xlab='CD133', main='Genes coexpressed with CD133')
# correlations = lm(dat[,'SOX11'] ~ dat[,'PROM1'])
# abline(correlations, col='red')
# summary(correlations)
# text(locator(1), 'R-squared: 0.1768')
# par(mfrow=c(1,1))


###########################################################################################################################


######################################### Build a heatmap of the correlated genes using the square matrix #####################

# Subset the dataframe with correlation values for those with high correlation and significance
prom1Cgenes = result[result[,2] > 0.1 & result[,4] < 0.05,]

# ALternative method that uses the standard deviation of the raw correlation values
corMean = mean(result[,1])
corSD = sd(result[,1])
prom1CgenesV2 = result[abs(result[,1]) > 3*corSD & result[,4] < 0.05,]

length(row.names(prom1Cgenes))
#write.table(prom1Cgenes, './manualCorrelation/140417_prom1CorrelatedGenes.txt', sep='\t')
length(row.names(prom1CgenesV2))
#write.table(prom1CgenesV2, './manualCorrelation/140417_prom1CorrelatedGenesUsed3SD.txt', sep='\t')

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
# write.table(squareAdjacency, '140430_squareAdjacencyMatrix.txt', sep='\t')

#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix,
#and calculate the corresponding dissimilarity:

# Calculate the topological overlap matrix
similarity = TOMsimilarity(squareAdjacency, TOMType='unsigned', verbose=3)
row.names(similarity) = row.names(squareAdjacency)
colnames(similarity) = row.names(squareAdjacency)
# write.table(similarity, '140430_similarSquareMatrix.txt', sep='\t')
###########################################################################################################################

######################################### Visulaising the network that was built ##########################################

# Calculate dissimilarity
dissTOM = 1-similarity

# Call the hierarchical clustering function. This is the faster implementation of WGCNA version
geneTree = flashClust(as.dist(dissTOM), method = "average")

#Branches of the dendrogram group together densely interconnected, highly co-expressed genes. Module identifcation amounts to the identification of individual branches
#There are several methods for branch cutting; our standard method is the Dynamic Tree Cut from the package dynamicTreeCut

# Module identification using dynamic tree cut. This is the most basic method and returns 3 modules when the cutHeight is 0.999 (default 0.99)
dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=0.999, method='tree')
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)


# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = plotTOM = dissTOM^6
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

# Plot the heatmap
TOMplot(plotTOM, geneTree, dynamicColors, main = "Network heatmap plot of CD133 coexpressed genes") #, terrainColors=FALSE)
        #,labRow=prom1CgenesNames, ColorsLeft=NA)

table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
###########################################################################################################################
