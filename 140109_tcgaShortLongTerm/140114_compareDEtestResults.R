setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/limmaResults/')
list.files()
# Read in the files into a list
affyFiles = list.files(pattern='*.affymetrix*')
agilentFiles = list.files(pattern='*.agilent*')
affy = lapply(affyFiles, read.delim, header=T)
names(affy) = affyFiles

agilFiles = lapply(agilentFiles, read.delim, header=T)
names(agilFiles) = agilentFiles

# Build a dataframe of the adjusted p-values from different methods
affyPval = cbind(affy[[1]]$adj.P.Val, affy[[2]]$adj.P.Val, affy[[3]]$adj.P.Val, affy[[4]]$adj.P.Val, affy[[5]]$adj.P.Val)
row.names(affyPval) = affy[[1]]$ID
colnames(affyPval) = affyFiles
affySig = -log10(affyPval)

# Build a dataframe of the adjusted p-values from different methods
agilentPval = cbind(agilFiles[[1]]$adj.P.Val, agilFiles[[2]]$adj.P.Val, agilFiles[[3]]$adj.P.Val)
row.names(agilentPval) = agilFiles[[1]]$ID
colnames(agilentPval) = agilentFiles
agilentSig = -log10(agilentPval)

head(affySig)
head(agilentSig)