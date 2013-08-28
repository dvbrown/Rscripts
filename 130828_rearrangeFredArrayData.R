#Arrange Fred's CSC microarray into a form amenable to GSEA
setwd('~/Documents/FredCSC/')
list.files(pattern='*.txt')

fdrPval = read.delim('Post-Hoc_PXR.txt',row.names=1)
avFC = read.delim('FC moyen PXR2+PXR6.txt', row.names=1)
logFC = read.delim('PXR-Transcriptome_Data.txt')