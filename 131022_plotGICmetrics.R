setwd('~/Documents/RNAdata/')

data = read.delim('RNAseqProgress.txt')
da = data[c(1:7),c(4,5,6,7,8,9,11,12,14)]
da[4,3] = as.factor('primary')

colors = c('aquamarine', 'cyan', 'skyblue', 'seagreen', 'palegreen', 'springgreen')

barplot(da$Survival..months., main='Survival', ylab='Survival (months)', names.arg=da$Clone, xlab='Patient ID')