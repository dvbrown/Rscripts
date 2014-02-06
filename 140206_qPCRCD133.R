# Analyse the the clones I was sequencing in batch 1 and compare primary and recurrent
source('~/Documents/Rscripts/qPCRFunctions.R')
################################## IO ################################## 
setwd('~/Documents/RNAdata/qPCRexpt/140206_cd133posNeg/')

# Convert the 384 well plate into the sample label
# Read in qPCR data
plate = transposeLinear('140206_384wellMap.txt')
# Split the transposed file into source and gene
plateMap = splitSampleName(plate)

cp = read.delim('140206_Cps.txt', skip=1)
tm = read.delim('140206_melt.txt', skip=1)

########################################################################

data = buildDataFrameForddCT(plateMap, cp)

# Return a datframe with the replicate data alongside
replicates = extractReplicates(c(1:384), data)
rawData = replicates[[3]]

rep1 = replicates[[1]]
rep2 = replicates[[2]]

# ################################### Munging the Tm manually ####################################
# tm = tm[,c(3,4,5,6)]
# mySampleLabels = sampleLabels[c(313:340),]
# tm = merge(mySampleLabels, tm, by.x='V1', by.y='Pos')
# #myTm = tm[c(313:340),]
# colnames(tm) = c('well', 'gene','sample', 'Tm1', 'Tm2' )
# repTm = extractReplicates(c(1:28), tm)
# repTm = repTm[[3]]
# Plot Tm
#plot(repTm$Tm1.x, repTm$Tm1.y, main='Replicate accuracy Tm', ylab='Tm')
#abline(lm(repTm$Tm1.x ~ repTm$Tm1.y), col='red')

# Set some graphs
par(las=2, mfrow=c(2,1))

# Plot the frequency of the Cp values
hist(data$Cp, 25, main='qPCR results in raw form', col='blue', xlab='Crossing point')

# Plot correlation of the replicates

plot(rawData$Cp.x, rawData$Cp.y, main='Replicate accuracy Cp', ylab='Cp replicate 1', xlab='Cp replicate 2')
abline(lm(rawData$Cp.x ~ rawData$Cp.y), col='red')
summary(lm(rawData$Cp.x ~ rawData$Cp.y))
text(50, 50, labels='R squared = 0.979')
