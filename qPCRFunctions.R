library(reshape)

# Intialise the package at the end by building a list containing all the functions in this script
#package.skeleton(name = 'qPCRcustomFunctions', list=c(), path='~/Documents/Rscripts/', force=T)

# Transpose the 384 well map (tab text) using a python script and output into another file
transposeLinear = function(well384Map, linearMapFile='output.txt') {
  pythonCall = paste('~/Documents/Eclipseworkspace/Bioinformatics/Filtering/transposeLinear.py -i', well384Map, '>', linearMapFile, sep=' ')
  system(pythonCall)
  sampleLabels = read.delim(linearMapFile, header=F)
  colnames(sampleLabels) = c('location', 'sample')
  return (sampleLabels)
}

#split sample names by whitespace. Takes a dataFrame and splits anything with whitespace into 2 columns
splitSampleName = function(plateMap) {
  # The column with sample is the vector containing the sample names you wish to split
  splitted = colsplit.factor(plateMap[['sample']], split = " ", names = c('origin', 'gene'))
  result = cbind(plateMap, splitted)
  return (result)
}

buildDataFrameForddCT = function(plateMap, CtData) {
  #Bind the dataframes containing the sample labels with the raw data itself and remove useless columns
  rawData = CtData[,c(3,4,5,8)]
  result = merge.data.frame(plateMap, rawData, by.x='location', by.y='Pos')
  return (result)
}

extractReplicates <- function (indexes, CtData) {
  indexes = c(1:384)
  
  # NEED TO CHECK THE NA CASES AS THE MERGE FUNCTION MAY DROP THEM!
  CtData = data#[(data[,3]) != <NA>,]
  # Subset each Cp into its replicates. Takes a vector with the indexes to to subset and then takes the
  # even entries and odd entries separately from the dataframe containing cp values
  even = indexes[indexes%%2 == 0]
  odd = indexes[indexes%%2 == 1]
  
  rep1 = CtData[odd, c(1:6)]
  rep1 = na.omit(rep1)
  rep2 = CtData[even, c(1:6)]
  rep2 = na.omit(rep2)
  boundData = merge(rep1, rep2, by.x='sample', by.y='sample')
  ################
  boundData = boundData[,c(1,2,3,4,6,7,11)]
  boundData$mean = rowMeans(cbind(boundData$Cp.x, boundData$Cp.y))
  result = list(rep1, rep2, boundData)
  return (result)
}