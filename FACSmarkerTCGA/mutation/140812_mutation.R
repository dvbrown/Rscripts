library(jsonlite)

setwd('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_mutation-2014-05-02/')
rawData = read.delim('genomicMatrix', row.names=1)
rawData[1,1] = 'blank'

row.names(rawData) = rawData$sample