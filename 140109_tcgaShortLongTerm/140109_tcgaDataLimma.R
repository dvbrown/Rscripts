# This script performs the differential expression testing
library(affy)
library(limma)
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
list.files()

overlapClinicalExpression = function(designMatrix, expresionDataFrame) {
  dSample = design[,1]
  eSample = colnames(expresionDataFrame)
  overlap = intersect(dSample, eSample)
  row.names(designMatrix) = designMatrix[,1]
  dOutput = designMatrix[overlap,]
  #dOutput = subset.data.frame(designMatrix, designMatrix[,1] == overlap)
  sOutput = expresionDataFrame[,overlap]
  output = list(dOutput, sOutput)
  return (output)
}

affy = read.delim('140109_affyMetrix.txt')

agilent = read.delim('140109_agilent.txt')
par(mfrow=c(2,1))
boxplot(affy[,c(1,5,10,50,100,150,200,211,333,444,500)], par(las=2))
boxplot(agilent[,c(1,5,10,50,100,150,200,211,333,444,499,501,512)], par(las=2))
design = readTargets('140109_targets.txt', sep='\t')

affyDesign = overlapClinicalExpression(design, affy)

# Affy was subjected gene_rma__data
affyEst = ExpressionSet(assayData=as.matrix(affyDesign[[2]]))

# Agilent data was => unc_lowess_normalization_gene_level__data
agilentEst = ExpressionSet(assayData=as.matrix(agilent))