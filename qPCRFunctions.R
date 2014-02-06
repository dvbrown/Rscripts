# Intialise the package at the end by building a list containing all the functions in this script

#package.skeleton(name = 'qPCRcustomFunctions', list=c(), path='~/Documents/Rscripts/', force=T)

# Transpose the 384 well map (tab text) using a python script and output into another file
transposeLinear = function(well384Map, linearMapFile='output.txt') {
  pythonCall = paste('~/Documents/Eclipseworkspace/Bioinformatics/Filtering/transposeLinear.py -i', well384Map, '>', linearMapFile, sep=' ')
  system(pythonCall)
  sampleLabels = read.delim(linearMapFile, header=F)
  return (sampleLabels)
}

#split sample names by whitespace. Takes a dataFrame and splits anything with whitespace into 2 columns
splitSampleName = function(aggregatedSample) {
  df = transform(aggregatedSample, sample = colsplit(sample, split = " ", names = c('sample', 'gene')))
  return (df)
}