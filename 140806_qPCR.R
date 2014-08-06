# Big qPCR expt

setwd('~/Documents/RNAdata/qPCRexpt/140728_mixedPop/140805_honourMashup/')
rawData = read.delim('140805_honoursMashup.txt', skip=1)
map = '../140805_384map.txt'

# Convert the 384 well plate into the sample label
plate = transposeLinear(map)
plateMap = splitSampleName(plate)
plateMap = plateMap[,c(1:4)]

data = buildDataFrameFromddCT(plateMap, rawData)
data = sort.dataframe(data, 1, highFirst=F)