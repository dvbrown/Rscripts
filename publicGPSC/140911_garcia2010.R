library(limma)
library(affy)
library(gplots)
library(ggplot2)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

##################### IO ###########################

setwd('~/Documents/public-datasets/GPSC_subgroups/garcia2010/GSE18015_RAW/')
list.files()

rawData = ReadAffy()
dm = readTargets('../designMatrix.txt')
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

par(mfrow=(c(2,1)))
# boxData = exprs(rawData)
# colnames(boxData) = dm$cellLine
boxplot(rawData, col=rainbow(21), main='Bhat unNormalised data', las=1, mar=c(5,5,5,2))
data = rma(rawData, verbose=T)
norm = exprs(data)
colnames(norm) = dm$cellLine
boxplot(norm, col=rainbow(21), main='Bhat normalised data', las=2, mar=c(5,5,5,2))
par(mfrow=(c(1,1)))