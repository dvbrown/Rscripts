# Analyse the the clones I was sequencing in batch 1 and compare primary and recurrent
source('~/Documents/Rscripts/qPCRFunctions.R')
################################## IO ################################## 
setwd('~/Documents/RNAdata/qPCRexpt/140206_cd133posNeg/')

# Convert the 384 well plate into the sample label
# Read in qPCR data
plate = transposeLinear('140206_384wellMap.txt')
# Split the transposed file into source and gene
plateMap = splitSampleName(plate)

cp = read.delim('140123_Pmd_Dvb.txt', skip=1)
tm = read.delim('140123_Pmd_DvbTmCALL.txt', skip=1)


########################################################################

buildCpDataframe <- function (sampleLabel, cp) {
  # Take the full dataframe listing the well and sample name and bind it to the file with well position and Cp
  mySampleLabels = sampleLabels[c(313:340),]
  myCp = cp[c(313:340),]
  # Merge labels with data
  data = merge(mySampleLabels, myCp, by.x='V1', by.y='Pos')
  data = data[,c(1,2,5,6,9)]
  colnames(data) = c('well', 'gene','sample', 'Cp', 'notes' )
  return (data)
}

extractReplicates <- function (indexes, cpData) {
  # Subset each Cp into its replicates. Takes a vector with the indexes to to subset and then takes the
  # even entries and odd entries separately from the dataframe containing cp values
  even = indexes[indexes%%2 == 0]
  odd = indexes[indexes%%2 == 1]
  rep1 = cpData[odd, ]
  rep2 = cpData[even, ]
  boundData = merge(rep1, rep2, by.x='gene', by.y='gene')
  boundData = boundData[,c(1,2,4,6,8)]
  boundData$mean = rowMeans(cbind(boundData$Cp.x, boundData$Cp.y))
  result = list(rep1, rep2, boundData)
  return (result)
}

################################### Munging the Tm manually ####################################
tm = tm[,c(3,4,5,6)]
mySampleLabels = sampleLabels[c(313:340),]
tm = merge(mySampleLabels, tm, by.x='V1', by.y='Pos')
#myTm = tm[c(313:340),]
colnames(tm) = c('well', 'gene','sample', 'Tm1', 'Tm2' )
repTm = extractReplicates(c(1:28), tm)
repTm = repTm[[3]]

# Set some graphs
par(las=2, mfrow=c(2,2))

# Plot Tm
plot(repTm$Tm1.x, repTm$Tm1.y, main='Replicate accuracy Tm', ylab='Tm')
abline(lm(repTm$Tm1.x ~ repTm$Tm1.y), col='red')

# Plot Cp
data = buildCpDataframe(sampleLabels, cp)
replicates = extractReplicates(c(1:28), data)
rep1 = replicates[[1]]
rep2 = replicates[[2]]
completeData = replicates[[3]]
completeData$mean = rowMeans(cbind(completeData$Cp.x, completeData$Cp.y))

plot(completeData$Cp.x, completeData$Cp.y, main='Replicate accuracy Cp', ylab='Cp')
abline(lm(completeData$Cp.x ~ completeData$Cp.y), col='red')
summary(lm(completeData$Cp.x ~ completeData$Cp.y))
text(50, 50, labels='R squared = 0.6964')

# barplot the mean of the Cp
barplot(completeData$mean, space=0.2, names.arg=completeData$gene, main='Mean of the Cp', col=rainbow(12),
        ylab='Crossing point')
abline(a=20, b=0)