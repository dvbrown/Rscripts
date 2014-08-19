source('~/Documents/Rscripts/multiplot.R')
source('~/Documents/Rscripts/qPCRFunctions.R')
setwd('~/Documents/RNAdata/qPCRexpt/140819_mixedPopAgain/')

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
write.table(rawData, '140403_rawData.txt', sep='\t')