library(ggplot2)
#this analyzes replicate data (eg qPCR)

raw.data = read.delim('yourfile')

na.data = raw.data$Cp
na.data = replace(na.data, is.na(na.data), '45') #this turns NAs into 40 Ct values
data = cbind(raw.data, na.data) #check the columns line up
data = as.data.frame(data)
data$na.data = as.numeric(data$na.data)

#calculate column mean and standard deviation
data.mean = aggregate(data$na.data, by=list(data$Sample), mean)
data.sd = aggregate(data$na.data, by=list(data$Sample), sd)
final.data = merge(data.mean, data.sd, by.x='Group.1', by.y='Group.1')

write.table(final.data,'date', sep='\t')

#select data with Ct below 40
good.data = subset.data.frame(final.data, 'mean.Ct' <40, select=c(Sample, mean.Ct, sd.Ct))

#encapsulate in a function
aggregateReplicates = function(dataFrame, dataColumn, sampleName) {
  data.mean = aggregate(dataColumn, by=list(sampleName), mean)
  data.sd = aggregate(dataColumn, by=list(sampleName), sd)
  final.data = merge(data.mean, data.sd, by.x='Group.1', by.y='Group.1')
  return (final.data)
}