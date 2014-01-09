# This script will import some TCGA data and test for differential expression between long and short term survivors
source('~/Documents/Rscripts/120704-sortDataFrame.R')
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/')
list.files()

makeDesignMatrix <- function (dataFrame,colNumSurvival, colNumAge, partitionSurvivalTime) {
  # Takes a dataframe containing survival data and patient names (the column names). Also specify which column contains survival data
  # Last argument is integer defining the cutoff (days) for short-term survivors (inclusive) and long-term survivors (exclusive)
  # Not sure how the row.names first line part is working. Check this if things go wrong
  patientNames = row.names(dataFrame)
  survival = dataFrame[,colNumSurvival]
  age = survival = dataFrame[,colNumAge]
  status = ifelse(survival, survival <= partitionSurvivalTime, survival > partitionSurvivalTime )
  status = (ifelse(status, 'short', 'long'))
  design = cbind(patientNames, age, status)
  design = as.data.frame(design)
  design$status = as.factor(design$status)
  return (design)
}

# This is normalised level 3 data from the TCGA. Downloaded using firehose and then filtered using python scripts

clinical = read.delim('140108_clinicalDataPart1.txt')
clinical2 = read.delim('140108_clinicalDataPart2.txt')
# Transpose clinical data (both parts are identical) to make it look like gene expression. Convert to correct type
clinical2 = t(clinical2)
row.names(clinical2) = colnames(clinical)
write.table(clinical2, '140109_clinicalDataTCGA.txt', sep='\t')

affy = read.delim('140108_affymetrixGem.txt', row.names=1)
agilent= read.delim('140108_agilentPart1gem.txt')
agilent2 = read.delim('140108_agilentPart2gem.txt')
# Agilent2 is out of order
agilent2 = sort.dataframe(agilent2, 1, highFirst=FALSE)

# Bind the 2 parts of agilent together
agilentTotal = merge.data.frame(agilent, agilent2, by.x='Hybridizatio', by.y='Hybridizatio')
row.names(agilentTotal) = agilentTotal[,1]
agilentTotal = agilentTotal[,c(2:513)]
rm(agilent, agilent2)

# Build the targets file (design matrix for use by limma)
# The vector containing clinical patient names

# Make a vector of short or long term survivors based on survival time
clinical2 = read.delim('140109_clinicalDataTCGA.txt')
# Check the row.names of clinical 2 for the below function to work. Should be patient names
design = makeDesignMatrix(clinical2, 1, 4,365)
write.table(design, '140109_targets.txt', sep='\t')