# Analyse the the clones I was sequencing in batch 1 and compare primary and recurrent
library(lattice)
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
row.names(rawData) = rawData$sample

rep1 = replicates[[1]]
rep2 = replicates[[2]]

par(las=2)
barchart(meanCP~gene.x,data=rawData,groups=origin.x, 
        scales=list(x=list(rot=90,cex=0.8)), main='Raw Cp values GIC cells')

# Set some graphs
par(las=2, mfrow=c(1,2))

# Plot the frequency of the Cp values
hist(data$Cp, 25, main='qPCR results in raw form', col='blue', xlab='Crossing point')

# Plot correlation of the replicates

plot(rawData$Cp.x, rawData$Cp.y, main='Replicate accuracy Cp', ylab='Cp replicate 1', xlab='Cp replicate 2')
abline(lm(rawData$Cp.x ~ rawData$Cp.y), col='red')
summary(lm(rawData$Cp.x ~ rawData$Cp.y))
text(50, 50, labels='R squared = 0.979')

# Generate the ddCt values by calling a custom function
rawData$ddCt_GAP_020_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                             houseKeepingGene='GAPDH', referenceSample='020_P', data=rawData)

rawData$ddCt_B2M_020_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                                       houseKeepingGene='B2M', referenceSample='020_P', data=rawData)

rawData$ddCt_GAP_030_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='030_P', data=rawData)
(rawData)


write.table(rawData, '140207_ddCt_calculations.txt', sep='\t', row.names=F)

# Subset the dataFrame by the comparisons of interest to make graphs clearer
shortLongSurvival = rawData[rawData$origin.x %in% c('020_P', '041_P'),]
primaryRecurrent = rawData[rawData$origin.x %in% c('030_P', '030a_P'),]
cd133negPos = rawData[rawData$origin.x %in% c('020_N','020_P', '030a_N', '030a_P', '041_N', '041_P'),]

# Generate ddCT values for the cd133 negative cell lines
cd133negPos$ddCt_020_N = ddCTcalculate(geneOfInterest=cd133negPos$gene.x, sampleOfInterest=cd133negPos$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='020_N', data=cd133negPos)
cd133negPos$ddCt_030_N = ddCTcalculate(geneOfInterest=cd133negPos$gene.x, sampleOfInterest=cd133negPos$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='030a_N', data=cd133negPos)
cd133negPos$ddCt_041_N = ddCTcalculate(geneOfInterest=cd133negPos$gene.x, sampleOfInterest=cd133negPos$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='041_N', data=cd133negPos)

# ######################################## Plot the ddCt values ####################################
# First compare B2M and GAPDH as house keeping genes
p1 = barchart(ddCt_B2M_020_P~origin.x,data=rawData,groups=gene.x, 
         scales=list(x=list(rot=90,cex=0.8)), main='B2M as the house keeping gene')
p2 = barchart(ddCt_GAP_020_P~origin.x,data=rawData,groups=gene.x, 
         scales=list(x=list(rot=90,cex=0.8)), main='GAPDH as the house keeping gene')

p3 = plot_ddCt(ddCt_B2M_020_P ~ origin.x, rawData,'Look no hands')

print(p1, position=c(0, .6, 1, 1), more=TRUE)
print(p2, position=c(0, 0, 1, .4))
#################################### Now the actual plots of interest
barchart(ddCt_B2M_020_P~origin.x,data=shortLongSurvival,groups=gene.x, 
         scales=list(x=list(rot=90,cex=0.8)), main='Short term vs long-term survivors')

s2 = plot_ddCt(ddCt_GAP_030_P~origin.x, primaryRecurrent, 'Primary vs recurrent tumours')
print(sl, position=c(0, .6, 1, 1), more=TRUE)
print(s2, position=c(0, 0, 1, .4))

s3 = plot_ddCt(ddCt_020_N~origin.x, dataFrame=cd133negPos,
               title='CD133 negative vs CD133 positive tumours')
print(s3)

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