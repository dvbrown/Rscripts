library(ggplot2)
#this analyzes replicate data (eg qPCR)

raw.data = read.delim('yourfile') 

analyseqPCR = function(dataFile) {
  data$Cp = replace(data$Cp, is.na(data$Cp), 45) #this turns NAs into 45 Ct values

  #calculate column mean and standard deviation
  data.mean = aggregate(data$Cp, by=list(data$sample), mean)
  data.sd = aggregate(data$Cp, by=list(data$sample), sd)
  final.data = merge(data.mean, data.sd, by.x='Group.1', by.y='Group.1')
  colnames(final.data) = c('sample', 'mean', 'stdDev')
  final.data  
  }
#select data with Ct below 40
good.data = subset.data.frame(final.data, mean <40, select=c(sample, mean, stdDev))