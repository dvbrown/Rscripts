
setwd('~/Documents/RNAdata/qPCRexpt/140331_cd133s/')

########################################## IO #######################################################

# Convert the 384 well plate into the sample label
# Read in qPCR data. Stupid python script has too many dependies and the files are everywhere
plate = system('../exec/transposeLinear.py -i ../140331_cd133s/384well_plateMap.txt > ../140331_cd133s/linear.txt')
plate = read.delim('140206_linearSamples.txt', header=F)
colnames(plate) = c('well', 'sample')
# Split the transposed file into source and gene
plateMap = splitSampleName(plate)

cp = read.delim('140403_cpScores.txt', skip=1)
tm = read.delim('140403_melt.txt', skip=1)

######################################### Filter and prepare data #####################################