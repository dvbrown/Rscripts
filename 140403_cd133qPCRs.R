source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/qPCRFunctions.R')
setwd('~/Documents/RNAdata/qPCRexpt/140331_cd133s/')

########################################## IO #######################################################

# Convert the 384 well plate into the sample label
# Read in qPCR data. Stupid python script has too many dependies and the files are everywhere
filename = '384well_plateMap.txt'
# Convert the 384 well plate into the sample label
# Read in qPCR data
plate = transposeLinear(filename)
# Split the transposed file into source and gene
plateMap = splitSampleName(plate)

cp = read.delim('140403_cpScores.txt', skip=1)
tm = read.delim('140403_melt.txt', skip=1)

######################################### Filter and prepare data #####################################