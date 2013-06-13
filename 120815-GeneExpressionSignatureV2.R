#A script to calculate similarity scores for gene expression
setwd('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/UNC__AgilentG4502A_07_1/')
source('~/Documents/Rscripts/120518-initialize.R')

geneMatrix = read.delim('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/UNC__AgilentG4502A_07_1/120814-CREBexpressionMatrixSorted.txt')
row.names(geneMatrix) = geneMatrix$genes
geneMatrix = geneMatrix[,2:92]

#Get the summary stats for the matrix
stats = summary(geneMatrix)

#In excel get the upper and lower quartiles and put into a vector for the CREB signature. Order alphabetically
sig = read.delim('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/UNC__AgilentG4502A_07_1/120813-CREBsignatureSorted.txt', header=FALSE)
signature = as.matrix(sig$V2)
row.names(signature) = sig$V1

#generate the distance matrix. Make sure that the order of genes in the expression matrix and signature are the same!
similarityMatrix = function(geneExpressionMatrix, GeneSignature) {
  index = 1
  distanceMatrix = as.matrix(geneExpressionMatrix)
  signature = as.vector(GeneSignature)
  result = vector(mode='numeric', length = length(geneExpressionMatrix))
  while (index <= length(geneExpressionMatrix)) {
    result[index] = sqrt((sum(distanceMatrix[,index] - signature) ^ 2))
    index = index + 1
  }
  result = as.matrix(result)
  row.names(result) = colnames(geneExpressionMatrix)
  return (result)
}
distanceMatrix = similarityMatrix(geneMatrix, signature)
#Obtain the column sum and mean for the distance matrix
hist(distanceMatrix, xlim=c(0, 20), ylim=c(0,20), breaks='fd', freq=T, xlab='Sum distance from CREB signature', ylab='Frequency', main='TCGA Agilent 1 sample similarity to 38 gene CREB signature', col='blue')

survival = read.delim('~/Documents/public-datasets/TCGA/clinicalData/120731-survivalDataStructered.txt')
survivalTime = survival[,c(4,5)]

#bind the survival data to the gene signature score data. This function needs debugging
bindSignatureSurvival = function(signatureScore, survivalTime) {
  distPaitents = row.names(signatureScore) #subset the survival data for paitents in the expression matrix
  sur = survivalTime[distPaitents,]
  signatureSurvival = merge(as.data.frame(signatureScore), sur, by.x = 'row.names', by.y='row.names')
  colnames(signatureSurvival) = c('Paitent', 'CREB.score', 'Survival', 'Censored')
  return (signatureSurvival)
}
distanceSurvival = bindSignatureSurvival(distanceMatrix, survivalTime)
#draw up some plots of the data
plot(distanceSurvival$Survival, distanceSurvival$CREB.score)