library(affyPLM)
library(limma)

setwd('~/Documents/public-datasets/rembrandt/rembrandt_GBM/')
allFiles = list.files()
allFiles 

# rawData = ReadAffy('00518392_U133P2.CEL', rm.mask=T, rm.outliers=T)
# exprData = exprs(rawData)
# rmaData = justRMA('00518392_U133P2.CEL')

rmaData = just(allFiles, rm.outliers=T, rm.mask=T, verbose=T, sampleNames=allFiles, destructive=T)