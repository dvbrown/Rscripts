library(reshape)
library(lattice)

# Intialise the package at the end by building a list containing all the functions in this script

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

extractReplicates <- function (indexes, ctData) {
  # Retreive the indexes of the 384 wellplate
  indexes = c(1:384)
  
  # Keep only cases with data in them as the merge function doesn't work with NAs
  CtData = data[complete.cases(data[,3]),]
  # Subset each Cp into its replicates. Takes a vector with the indexes to to subset and then takes the
  # even entries and odd entries separately from the dataframe containing cp values
  even = indexes[indexes%%2 == 0]
  odd = indexes[indexes%%2 == 1]
  
  rep1 = CtData[odd, c(1:6)]
  rep1 = na.omit(rep1)
  rep2 = CtData[even, c(1:6)]
  rep2 = na.omit(rep2)
  boundData = merge(rep1, rep2, by.x='sample', by.y='sample')
  ################ Remove columns that do not add information
  usefulData = boundData[,c(1,2,3,4,6,7,11)]
  # Compute the mean and the standard deviation of the replicates
  usefulData$meanCP = rowMeans(cbind(usefulData$Cp.x, usefulData$Cp.y))
  usefulData$stdDevCP = apply(cbind(usefulData$Cp.x, usefulData$Cp.y), 1, sd)
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

plot_ddCt = function(Expressionformula, dataFrame, title='A grouped barchart') {
  # This will make barcharts without error bars
  # Expression formula is of the type ddCt ~ cell type or whatever you want the bars to be grouped by
  p = barchart(Expressionformula, data = dataFrame, groups = gene.x, 
<<<<<<< HEAD
<<<<<<< HEAD
                scales = list(x = list(rot=90,cex=0.8)), main = title,  ylim=c(0,10),
               ylab=Expression)
=======
                scales = list(x = list(rot=90,cex=0.8)), main = title)
>>>>>>> parent of e21c6dc... the bar plot doesn't work anymore
=======
                scales = list(x = list(rot=90,cex=0.8)), main = title)
>>>>>>> parent of e21c6dc... the bar plot doesn't work anymore
  #returns a plot object that when you look at it plots stuff
  return (p)
}
# Run this at the end to intialise the package
#package.skeleton(name = 'qPCRcustomFunctions', list=c(), path='~/Documents/Rscripts/', force=F)