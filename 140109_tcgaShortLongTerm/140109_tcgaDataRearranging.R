# This script will import some TCGA data and test for differential expression between long and short term survivors
source('~/Documents/Rscripts/120704-sortDataFrame.R')
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/')
list.files()

clinical = read.delim('140108_clinicalDataPart1.txt')
clinical2 = read.delim('140108_clinicalDataPart2.txt', skip=0)
clinical2 = t(clinical2)

affy = read.delim('140108_affymetrixGem.txt', row.names=1)
agilent= read.delim('140108_agilentPart1gem.txt')
agilent2 = read.delim('140108_agilentPart2gem.txt')
# Agilent2 is out of order
agilent2 = sort.dataframe(agilent2, 1, highFirst=FALSE)

# Bind the 2 parts of agilent together
agilentTotal = merge.data.frame(agilent, agilent2, by.x='Hybridizatio', by.y='Hybridizatio')
row.names(agilentTotal) = agilentTotal[,1]
agilentTotal = agilentTotal[,c(2:513)]
rm(agilent, agilent2)