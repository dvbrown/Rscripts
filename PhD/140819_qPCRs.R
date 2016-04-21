setwd('~/Documents/RNAdata/qPCRexpt/140819_mixedPopAgain/')
source('~/Documents/Rscripts/multiplot.R')
source('~/Documents/Rscripts/qPCRFunctions.R')

extractReplicates <- function (indexes, ctData) {
    # Retreive the indexes of the 384 wellplate
    indexes = c(1:384)
    
    # Keep only cases with data in them as the merge function doesn't work with NAs
    CtData = ctData[complete.cases(ctData[,3]),]
    # Subset each Cp into its replicates. Takes a vector with the indexes to to subset and then takes the
    # even entries and odd entries separately from the dataframe containing cp values
    even = as.character(indexes[indexes%%2 == 0])
    odd = as.character(indexes[indexes%%2 == 1])
    even = paste('Sample', even, sep=' ')
    odd = paste('Sample', odd, sep=' ')
    
    #rep1 = CtData[odd, c(1:6)]
    rep1 = CtData[CtData$Name %in% odd, c(1:6)]
    
    #rep1 = rep1[complete.cases(rep1$sample),]
    
    #rep2 = CtData[even, c(1:6)]
    rep2 = CtData[CtData$Name %in% even, c(1:6)]
    
    #rep2 = rep2[complete.cases(rep2$sample),]
    
    boundData = merge(rep1, rep2, by.x='sample', by.y='sample')
    ################ Remove columns that do not add information
    usefulData = boundData[,c(1,2,3,4,6,7,11)]
    # Compute the mean and the standard deviation of the replicates
    usefulData$meanCP = rowMeans(cbind(usefulData$Cp.x, usefulData$Cp.y), na.rm=T)
    usefulData$stdDevCP = apply(cbind(usefulData$Cp.x, usefulData$Cp.y), 1, sd, na.rm=T)
    # Package the output in a list
    result = list(rep1, rep2, usefulData)
    return (result)
}

ddCTcalculate = function(geneOfInterest, sampleOfInterest='020_N', houseKeepingGene='GAPDH', referenceSample='020_N', data=rawData) {
    # The gene of interest and the sample of interest are both vectors from the dataFrame.
    # House keeping gene and 
    
    # This function has been successfully debugged and has been shown to work.
    sampleHouse = paste(sampleOfInterest, houseKeepingGene)
    SampleGene = paste(sampleOfInterest, geneOfInterest)
    # Extract the Cp of the hosue keeping gene
    houseCp = data[sampleHouse, 'meanCP']
    # place holder 020_N ATP5G3 as the current row
    geneCp = data[SampleGene, 'meanCP']
    # dCt calculation for the sample of interest
    dCt = houseCp - geneCp
    
    # Extract the meanCP for the reference sample. First get the index of the housekeeping gene, then the gene of interest
    refDctRowHouse = paste(referenceSample, houseKeepingGene)
    refDctRowGene = paste(referenceSample, geneOfInterest)
    # Calculate dCt for the reference sample
    referenceSample_dCt = data[refDctRowHouse, 'meanCP'] - data[refDctRowGene, 'meanCP']
    
    # Calculate ddCt
    ddCtNotSquared = dCt - referenceSample_dCt
    ddCt = 2^ddCtNotSquared
    return (ddCt)
}


source('~/Documents/Coding/workCode/Rscripts/qPCRFunctions.R')
source('~/Documents/Coding/workCode/Rscripts/multiplot.R')

# Convert the 384 well plate into the sample label
# Read in qPCR data. Stupid python script has too many dependies and the files are everywhere
filename = './140819_384map.txt'
# Convert the 384 well plate into the sample label
# Read in qPCR data
plate = transposeLinear(filename)
# Check the 384 well plate plain text file using cat or less. Gremlins in the file like empty wells and newlines will screw up the script

# Split the transposed file into source and gene
plateMap = splitSampleName(plate)

cp = read.delim('140819_cpScores.txt', skip=1)
tm = read.delim('140819_melt.txt', skip=1)

######################################### Filter and prepare data #####################################
data = buildDataFrameFromddCT(plateMap, cp)

# Return a datframe with the replicate data alongside
replicates = extractReplicates(c(1:384), data)
rawData = replicates[[3]]
row.names(rawData) = rawData$sample
write.table(rawData, '140819_mungedData.txt', sep='\t', row.names=F)

##### I started from home here ####
mungedData = read.delim('~/Dropbox//140819_mixedPopAgain//140819_rawData.txt', row.names=1)
mungedData = mungedData[c(1:35),]

# Plot correlation of the replicate
correlation = ggplot(data=mungedData, aes(x=Cp.y, y=Cp.x, color=gene.x)) + 
                geom_point(shape=19) + geom_smooth(method=lm, colour='red') +
                xlab("Replicate 1") + ylab("Replicate 2") + # Set axis labels
                ggtitle("Correlation of technical replicates") +  # Set title
                theme_bw(base_size=18)

############################################ Calculate the ddCt scores #################################################
# Subset the data by cell line
clones = levels(mungedData$origin.x)
clones = droplevels(as.factor(clones))

c035 = mungedData[mungedData$origin.x %in% c("035_DN", "035_44", "035_33", "035_DP"),]
mixed = mungedData[mungedData$origin.x %in% c("039_M", "020_M", "041_M"),]

# Change Cp values of 40 to NA
c035$meanCP[c035$meanCP == 40] = NA
mixed$meanCP[mixed$meanCP == 40] = NA

# The ddCt
c035$ddCt = ddCTcalculate(geneOfInterest=c035$gene.x, sampleOfInterest=c035$origin.x,
                          houseKeepingGene='GAPDH', referenceSample='035_DN', data=c035) 

mixed$ddCt = ddCTcalculate(geneOfInterest=mixed$gene.x, sampleOfInterest=mixed$origin.x,
                           houseKeepingGene='GAPDH', referenceSample='020_M', data=mixed)

mixedPop = ggplot(data=mixed, aes(x=gene.x, y=ddCt, fill=origin.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Sample") +      # Set legend title
    #scale_y_continuous(breaks = round(seq(min(bindData$ddCt), max(bindData$ddCt), by = 0.5),0.5)) + # This modifies the scale of the y axis.
    xlab("Gene") + ylab("Gene expression normalised\nto GPSC #035") + # Set axis labels
    ggtitle("Expression of Proneural/ Mesenchymal markers \nfor mixed population GPSC") +  # Set title
    theme_bw(base_size=14)
mixedPop + theme(axis.text.x = element_text(angle = 90, hjust = 1))

clone035 = ggplot(data=c035, aes(x=gene.x, y=ddCt, fill=origin.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Subpopulation") +      # Set legend title
    #scale_y_continuous(breaks = round(seq(min(bindData$ddCt), max(bindData$ddCt), by = 0.5),0.5)) + # This modifies the scale of the y axis.
    xlab("Gene") + ylab("Gene expression normalised \nto CD44-/CD133- subpopulation") + # Set axis labels
    ggtitle("Expression of Proneural/ Mesenchymal markers \nfor GPSC number #035") +  # Set title+
    theme_bw(base_size=14)
clone035 + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(clone035 + theme(axis.text.x = element_text(angle = 90, hjust = 1)), mixedPop + theme(axis.text.x = element_text(angle = 90, hjust = 1)),correlation, cols=2)