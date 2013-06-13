#Complete analysis using functions
setwd('~/Documents/RNAdata/qPCRexpt/')
require(reshape)
library(ggplot2)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

mapRawData = function(map, data) {
  #Take a vertical linear mapping file generated from the 384 well map with sample names and bind it to raw data from lightcycler
  map = t(map)
  map = as.vector(map)
  wellPos = read.delim('~/Documents/RNAdata/qPCRexpt/120817-384welllayouttranspose.txt', header=F)
  wellPos = as.character(wellPos[,1])
  wellMap = as.data.frame(cbind(wellPos[1:length(map)], map))
  annotatedData = merge.data.frame(wellMap, data, by.x='V1', by.y='Pos')
  row.names(annotatedData) = annotatedData[,1]
  return (annotatedData)
}
#plot raw data
plotRawCt = function(mappedData) {
  par(cex=0.5, las=2, cex.axis=0.8)
  thePlot = barplot(mappedData$Cp, space =1, names.arg=mappedData[,2], col=rainbow(length(mappedData$Cp)), main='Raw Ct values', ylab='Ct', xlab='sample', ylim=c(0,45))
  return (thePlot)
}
#get the replicates together
aggregateReplicates = function(dataFrame, dataColumn, sampleName) {
  #the dataColumn and sample name need to be explicit eg dataFrame$data
  data.mean = aggregate(dataColumn, by=list(sampleName), mean)
  data.sd = aggregate(dataColumn, by=list(sampleName), sd)
  final.data = merge(data.mean, data.sd, by.x='Group.1', by.y='Group.1')
  colnames(final.data) = c('sample', 'mean.Cp', 'sd.Cp')
  return (final.data)
}

#split sample names by whitespace
splitSampleName = function(aggregatedSample) {
  df = transform(aggregatedSample, sample = colsplit(sample, split = " ", names = c('sample', 'gene')))
  return (df)
}

#fix strange encoding of the splitSampleName function.
fixSplitSample = function(dataFrame, sample, Cp) {
  sample = dataFrame$sample
  x = cbind(sample, dataFrame$Cp, dataFrame$well)
  return(x)
}

#this function needs more debugging
errorBarPlot = function(dataFrame, sample, mean, sd) {
  thePlot=  ggplot(dataFrame) +
    geom_bar(aes(x=sample, y=mean) , position=position_dodge(width=1), fill='blue', colour='black') +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), 
                width=0.2,
                position=position_dodge(.9), colour='black') +
                  ylab('Mean Ct') + opts(axis.text.x = theme_text(angle=90, hjust=1.2, size=12), 
                                         title='qRT-PCR stem cells')
  return(thePlot)
}