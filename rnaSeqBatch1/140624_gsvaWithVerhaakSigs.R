# Try and call subtype with GSVA
library(GSVA)
library(biomaRt)

setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')
verhaakSig = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')
verhaakSigMung = list(verhaakSig[,1], verhaakSig[,2], verhaakSig[,3], verhaakSig[,4])
names(verhaakSigMung) = colnames(verhaakSig)

dataCorrected = as.matrix(read.delim('140213_normalisedCPM_facs.txt'))
dataLog= as.matrix(read.delim('140213_normalisedLog_CPM_facs.txt'))

# extract the gene names for the ensembl IDs

result = gsva(dataCorrected, verhaakSigMung)