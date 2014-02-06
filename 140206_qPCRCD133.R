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
noNAs = na.omit(data)
hist(noNAs$Cp, 25, main='qPCR results in raw form', col='yellow')

# Return a datframe with the replicate data alongside
replicates = extractReplicates(c(1:384), data)
data = replicates[[1]]

# ################################### Munging the Tm manually ####################################
# tm = tm[,c(3,4,5,6)]
# mySampleLabels = sampleLabels[c(313:340),]
# tm = merge(mySampleLabels, tm, by.x='V1', by.y='Pos')
# #myTm = tm[c(313:340),]
# colnames(tm) = c('well', 'gene','sample', 'Tm1', 'Tm2' )
# repTm = extractReplicates(c(1:28), tm)
# repTm = repTm[[3]]

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