setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/limmaResults/')
list.files()
# Read in the files into a list
affyFiles = list.files(pattern='*.affymetrix*')
agilentFiles = list.files(pattern='*.agilent*')
affy = lapply(affyFiles, read.delim, header=T)
names(affy) = affyFiles
agilFiles = lapply(agilentFiles, read.delim, header=T)
names(agilFiles) = agilentFiles