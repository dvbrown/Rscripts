# This script will take Agilent gene expression from the TCGA and attempt to subset it for survival

setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/140214_testSignaturesAgilent/')
list.files()

agilent = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/140110_agilentNoNulls.txt')
# The noNulls refers to being put through a python script where NULL was turned into NA

# Patient number 500 is weirdly high in agilent. Try this analysis with and witahout including it and potentially other samples
#agilent = agilent[,c(1:499,501:512)]
#boxplot(agilent[,c(1,5,10,50,100,150,200,211,333,444,499,500,501,511)], par(las=2, cex=0.8), main='Agilent Lowess normalised')

# Check that the data has been z-transformed by looking at a few genes
par(mfrow=c(2,2))
geneMatrix = t(agilent)
hist((geneMatrix[,'TP53']), main='p53 distribution', xlab='TP53')
hist((geneMatrix[,'GAPDH']), main='GAPDH distribution', xlab='GAPDH')
hist((geneMatrix[,'EGFR']), main='EGFR distribution', xlab='EGFR')
hist((geneMatrix[,'ACTB']), main='B-Actin distribution', xlab='B-Actin')
par(mfrow=c(1,1))
rm(geneMatrix) 