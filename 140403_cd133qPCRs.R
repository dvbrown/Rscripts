source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/qPCRFunctions.R')
setwd('~/Documents/RNAdata/qPCRexpt/140331_cd133s/')

subsetClones = function(cloneName) {
    # Takes the name of the clone as a string
    neg = paste(cloneName, '_neg', sep='')
    pos = paste(cloneName, '_pos', sep='')
    result = rawData[rawData$origin.x %in% c(neg, pos),]
    return (result)
}

########################################## IO #######################################################

# Convert the 384 well plate into the sample label
# Read in qPCR data. Stupid python script has too many dependies and the files are everywhere
filename = '384well_plateMap.txt'
# Convert the 384 well plate into the sample label
# Read in qPCR data
plate = transposeLinear(filename)
# Check the 384 well plate plain text file using cat or less. Gremlins in the file like empty wells and newlines will screw up the script

# Split the transposed file into source and gene
plateMap = splitSampleName(plate)

cp = read.delim('140403_cpScores.txt', skip=1)
tm = read.delim('140403_melt.txt', skip=1)

######################################### Filter and prepare data #####################################
data = buildDataFrameFromddCT(plateMap, cp)

# Return a datframe with the replicate data alongside
replicates = extractReplicates(c(1:384), data)
rawData = replicates[[3]]
row.names(rawData) = rawData$sample
write.table(rawData, '140403_rawData.txt', sep='\t')

# Extract replictes from the data set
rep1 = replicates[[1]]
rep2 = replicates[[2]]

par(las=2)
hist(data$Cp, main='Raw Cp values GIC cells', col='orange')

# Set some graphs
par(las=2, mfrow=c(1,2))

# Plot correlation of the replicates

plot(rawData$Cp.x, rawData$Cp.y, main='Replicate accuracy Cp', ylab='Cp replicate 1', xlab='Cp replicate 2', pch=16)
abline(lm(rawData$Cp.x ~ rawData$Cp.y), col='red')
summary(lm(rawData$Cp.x ~ rawData$Cp.y))
text(locator(1), labels='R squared = 0.9809')

############################################ Calculate the ddCt scores #################################################
# Subset the data by cell line
c011 = subsetClones('011')
c020 = subsetClones('020')
c041 = subsetClones('041')
c035 = subsetClones('035')
c030a = subsetClones('030a')
c020 = subsetClones('020')
c030 = subsetClones('030')

# The ddCt
c011$ddCt = ddCTcalculate(geneOfInterest=c011$gene.x, sampleOfInterest=c011$origin.x,
                                       houseKeepingGene='GAPDH', referenceSample='011_neg', data=c011)

c020$ddCt = ddCTcalculate(geneOfInterest=c020$gene.x, sampleOfInterest=c020$origin.x,
                     houseKeepingGene='GAPDH', referenceSample='020_neg', data=c020)

c030$ddCt = ddCTcalculate(geneOfInterest=c030$gene.x, sampleOfInterest=c030$origin.x,
                          houseKeepingGene='GAPDH', referenceSample='030_neg', data=c030)

c030a$ddCt = ddCTcalculate(geneOfInterest=c030a$gene.x, sampleOfInterest=c030a$origin.x,
                           houseKeepingGene='GAPDH', referenceSample='030a_neg', data=c030a)

c035$ddCt = ddCTcalculate(geneOfInterest=c035$gene.x, sampleOfInterest=c035$origin.x,
                          houseKeepingGene='GAPDH', referenceSample='035_neg', data=c035)

c041$ddCt = ddCTcalculate(geneOfInterest=c041$gene.x, sampleOfInterest=c041$origin.x,
                          houseKeepingGene='GAPDH', referenceSample='041_neg', data=c041)

# Bind the data together
bindData = rbind(c011, c030a, c035, c041, c020, c030)

######################################### Plot the ddCt values ########################################################
ls2 = ggplot(data=shortLongSurvival[!shortLongSurvival$gene.x %in% c('B2M','GAPDH','OCT4','HAPLN1','GFAP','ATP5G3'),], 
             aes(x=origin.x, y=ddCt_B2M_020_P, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("Gene expression normalised to long-term survivor") + # Set axis labels
    ggtitle("Follow up of genes identified as differentially expressed by RNA-seq") +  # Set title
    theme_bw(base_size=20)