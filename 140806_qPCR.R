# Big qPCR expt
source('~/Documents/Rscripts/qPCRFunctions.R')

setwd('~/Documents/RNAdata/qPCRexpt/140728_mixedPop/140805_honourMashup/')
# rawData = read.delim('140805_honoursMashup.txt', skip=1)
# map = '../140805_384map.txt'
# 
# # Convert the 384 well plate into the sample label
# plate = transposeLinear(map)
# plateMap = splitSampleName(plate)
# plateMap = plateMap[,c(1:4)]
# 
# data = buildDataFrameFromddCT(plateMap, rawData)
# Use the mixed order function to sort the dataframe by sample number
dataSorted = data[mixedorder(data$Name),]
#write.table(dataSorted, '140805_Cpvalues.txt', row.names=F, sep='\t')

dataMy = dataSorted[c(1:168),]
rm(data, dataSorted, plate, plateMap, rawData)

# Extract replicates
replicates = extractReplicates(c(1:168), dataMy)
data = replicates[[3]]

#write.table(data, '140808_replicatesDan.txt', sep='\t')
data = read.delim('output/140808_replicatesDanfixedNames.txt', row.names=1)
##################################### Some QC plots #####################################
plot(dataSorted$Cp, dataSorted$location, pch=16, col='pink', main='Raw Cp scores qRT-PCR expt with honours students')

# Plot correlation of the replicates
plot(data$Cp.x, data$Cp.y, main='Replicate accuracy Cp', ylab='Cp replicate 1', xlab='Cp replicate 2', pch=16)
abline(lm(data$Cp.x ~ data$Cp.y), col='red')
summary(lm(data$Cp.x ~ data$Cp.y))
text(locator(1), labels='R squared = 0.975')
par(mfrow=c(1,1))

############################################ Calculate the ddCt scores #################################################
# Subset the data by cell line
clones = levels(data$origin.x)
c035 = data[data$origin.x %in% c("035_33", "035_44", "035_DN", "035_DP", "035_M"),]
c041 = data[data$origin.x %in% c("041_3N", "041_3P", "041_M"),]
mixed = c041 = data[data$origin.x %in% c("020_M", "035_M","039_M","041_M" ),]

# The ddCt
c035$ddCt = ddCTcalculate(geneOfInterest=c035$gene.x, sampleOfInterest=c035$origin.x,
                          houseKeepingGene='GAPDH', referenceSample='035_M', data=c035)