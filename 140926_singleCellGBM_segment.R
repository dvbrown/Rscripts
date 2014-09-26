setwd('~/Documents/public-datasets/RNA-seq/anoop2014_singleCellGBM/')

sigData = read.delim('140926_signatureScoresAllPat.txt')
geneData = read.delim('GSE57872_GBM_data_matrix.txt', row.names=1)