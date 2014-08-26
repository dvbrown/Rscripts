library(limma)
library(affy)

setwd('~/Documents/public-datasets/GPSC_subgroups/bhat2013/GSE49009_RAW/')
list.files()

rawData = ReadAffy()
dm = readTargets('../designMatrix.txt')

par(mfrow=(c(2,1)))
boxplot(rawData, col=rainbow(21), main='Bhat unNormalised data', las=2, mar=c(5,5,5,2))
# data = rma(rawData, verbose=T)
# norm = exprs(data)
boxplot(norm, col=rainbow(21), main='Bhat normalised data', las=2, mar=c(5,5,5,2))
