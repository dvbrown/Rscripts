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

######################################### Filter and prepare data #####################################

data = buildDataFrameForddCT(plateMap, cp)

# Return a datframe with the replicate data alongside
replicates = extractReplicates(c(1:384), data)
rawData = replicates[[3]]
row.names(rawData) = rawData$sample

# Get rid of the primers that didn't work and remove their factor levels
rawData = rawData[rawData$gene.x %in% c('ATP5G3', 'B2M', 'CREB1', 'GAPDH', 'GFAP', 'GRIA2', 'HAPLN1',
                                        'LAMB1', 'OCT4', 'PROM1', 'REST'),]
rawData = droplevels(rawData)

# Extract replictes from the data set
rep1 = replicates[[1]]
rep2 = replicates[[2]]

# par(las=2)
# barchart(meanCP~gene.x,data=rawData,groups=origin.x, 
#         scales=list(x=list(rot=90,cex=0.8)), main='Raw Cp values GIC cells')
# 
# # Set some graphs
# par(las=2, mfrow=c(1,2))
# 
# # Plot the frequency of the Cp values
# hist(data$Cp, 25, main='qPCR results in raw form', col='blue', xlab='Crossing point')
# 
# # Plot correlation of the replicates
# 
# plot(rawData$Cp.x, rawData$Cp.y, main='Replicate accuracy Cp', ylab='Cp replicate 1', xlab='Cp replicate 2')
# abline(lm(rawData$Cp.x ~ rawData$Cp.y), col='red')
# summary(lm(rawData$Cp.x ~ rawData$Cp.y))
# text(50, 50, labels='R squared = 0.979')

##########################################  Generate the ddCt values by calling a custom function ######################################### 
rawData$ddCt_GAP_020_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                             houseKeepingGene='GAPDH', referenceSample='020_P', data=rawData)

rawData$ddCt_B2M_020_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                                       houseKeepingGene='B2M', referenceSample='020_P', data=rawData)

rawData$ddCt_GAP_030_P = ddCTcalculate(geneOfInterest=rawData$gene.x, sampleOfInterest=rawData$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='030_P', data=rawData)
(rawData)


#write.table(rawData, '140207_ddCt_calculations.txt', sep='\t', row.names=F)
######################################### ######################################### ######################################### 

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
system('mkdir qcplots')

pdf("qcplots/B2M.pdf")
barchart(log2(ddCt_B2M_020_P)~origin.x,data=rawData,groups=gene.x, 
         scales=list(x=list(rot=90,cex=0.8)), main='B2M as the house keeping gene',
              auto.key=list(space="top", columns=4,
                            title="gene legend", cex.title=1))
dev.off()

pdf('140210_GAPDH_housekeeper.pdf')
barchart(log2(ddCt_GAP_020_P)~origin.x,data=rawData,groups=gene.x, 
         scales=list(x=list(rot=90,cex=0.8)), main='GAPDH as the house keeping gene',
              auto.key=list(space="top", columns=4,
                            title="genes", cex.title=1))
dev.off()


# Use custom function from now on
p3 = plot_ddCt(ddCt_B2M_020_P ~ origin.x, rawData,'Look no hands')

print(p1, position=c(0, .6, 1, 1), more=TRUE)
print(p2, position=c(0, 0, 1, .4))
#################################### Now the actual plots of interest
lSplot = barchart((ddCt_B2M_020_P)~origin.x,data=shortLongSurvival,groups=gene.x, ylab='ddCt',
         scales=list(x=list(rot=90,cex=0.8)), main='Short term vs long-term survivors', 
         auto.key=list(space="top", columns=4,
                       title="genes", cex.title=1))
update(lSplot, par.settings = list(fontsize = list(text = 18, points = 4)))


# Changed to logs
sl = plot_ddCt(log2(ddCt_B2M_020_P)~origin.x, shortLongSurvival, 'Short vs Long term survival')

s2 = plot_ddCt(log2(ddCt_GAP_030_P)~origin.x, primaryRecurrent, 'Primary vs recurrent tumours', yaxisLabel='ddCT')
update(s2, par.settings = list(fontsize = list(text = 18, points = 4)))

print(sl, position=c(0, .6, 1, 1), more=TRUE)
print(s2, position=c(0, 0, 1, .4))

######################### Build a dataframe for individual CD133 ddCT then rbind for plotting ########
# 030a
cd133_30a = rawData[c(33:53),c(1:9)]
cd133_30a$ddCt = ddCTcalculate(geneOfInterest=cd133_30a$gene.x, sampleOfInterest=cd133_30a$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='030a_N', data=cd133_30a)

#041
cd133_41 = rawData[c(54:75),c(1:9)]
cd133_41$ddCt = ddCTcalculate(geneOfInterest=cd133_41$gene.x, sampleOfInterest=cd133_41$origin.x,
                                     houseKeepingGene='GAPDH', referenceSample='041_N', data=cd133_41)

#020
cd133_20 = rawData[c(1:22),c(1:9)]
cd133_20$ddCt = ddCTcalculate(geneOfInterest=cd133_20$gene.x, sampleOfInterest=cd133_20$origin.x,
                              houseKeepingGene='GAPDH', referenceSample='020_N', data=cd133_20)

cd133negPos = rbind(cd133_20, cd133_30a, cd133_41)

cd133Plot = barchart(log(ddCt)~origin.x,data=cd133negPos,groups=gene.x, ylab='ddCt',
                  scales=list(x=list(rot=90,cex=0.8)), main='CD133 negative vs positive', 
                  auto.key=list(space="top", columns=3,
                                title="genes", cex.title=1))
update(cd133Plot, par.settings = list(fontsize = list(text = 18, points = 4)))


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