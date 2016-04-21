#A script to calculate similarity scores for gene expression
library(proxy)
setwd('Documents/public-datasets/TCGA/expressionData/Expression-Genes/UNC__AgilentG4502A_07_1/')
source('~/Documents/Rscripts/120518-initialize.R')

geneMatrix = read.delim('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/UNC__AgilentG4502A_07_1/120814-CREBexpressionMatrixSorted.txt')
row.names(geneMatrix) = geneMatrix$genes
geneMatrix = geneMatrix[,2:92]

#Get the summary stats for the matrix
stats = summary(geneMatrix)

#In excel get the upper and lower quartiles and put into a vector for the CREB signature. Order alphabetically
sig = read.delim('120813-CREBsignatureSorted.txt', header=FALSE)
signature = as.matrix(sig$V2)
row.names(signature) = sig$V1

#generate the distance matrix. Make sure that the order of genes in the expression matrix and signature are the same!
distanceMatrix = function(geneExpressionMatrix, GeneSignature) {
  index = 1
  distanceMatrix = geneExpressionMatrix[,,]
  while (index <= length(geneExpressionMatrix)) {
    distanceMatrix[,index] = dist(geneExpressionMatrix[,index], GeneSignature, by_rows=FALSE, pairwise=TRUE) #distance measure
    index = index + 1
  }
  return (distanceMatrix)
}

similarityMatrix = function(geneExpressionMatrix, GeneSignature) {
  index = 1
  distanceMatrix = geneExpressionMatrix[,,]
  while (index <= length(geneExpressionMatrix)) {
    distanceMatrix[,index] = simil(geneExpressionMatrix[,index], GeneSignature, pairwise=TRUE, by_rows=FALSE) #similarity measure
    index = index + 1
  }
  return (distanceMatrix)
}

#Obtain the column sum and mean for the distance matrix
meanDistance = colMeans(distanceMatrix)
sumDistance = colSums(distanceMatrix)
hist(sumDistance, xlim=c(10, 40), breaks='Scott', freq=T, xlab='Sum distance from CREB signature', ylab='Frequency', main='TCGA Agilent 1 sample similarity to 38 gene CREB signature', col='blue')
distanceSummary = summary(sumDistance)

#bind the survival data to the colSum similarity data
survival = read.delim('~/Documents/public-datasets/TCGA/clinicalData/120731-survivalDataStructered.txt')
sur = survival[,c(4,5)]
distPaitents = colnames(distanceMatrix)
sur1 = sur[distPaitents,]