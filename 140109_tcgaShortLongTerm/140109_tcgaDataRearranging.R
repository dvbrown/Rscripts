# This script will import some TCGA data and test for differential expression between long and short term survivors

setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/')
list.files()

clinical = read.delim('140108_clinicalDataPart1.txt')
clinical2 = read.delim('140108_clinicalDataPart2.txt', skip=0)
clinical2 = t(clinical2)

affy = read.delim('140108_affymetrixGem.txt')
agilent= read.delim('140108_agilentPart1gem.txt')
agilent2 = read.delim('140108_agilentPart2gem.txt')