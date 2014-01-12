# This script performs the differential expression testing
library(affy)
library(limma)
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
list.files()

overlapClinicalExpression = function(designMatrix, expresionDataFrame) {
  # Return a list with the first element is the overlapped design matrix. The second element is the gene expression matrix
  row.names(designMatrix) = designMatrix[,1]
  dSample = design[,1]
  eSample = colnames(expresionDataFrame)
  overlap = intersect(eSample, dSample)
  dOutput = designMatrix[overlap,]
  sOutput = expresionDataFrame[,overlap]
  output = list(dOutput, sOutput)
  return (output)
}

affy = read.delim('140110_affyNoNulls.txt')

agilent = read.delim('140110_agilentNoNulls.txt')
# Patient number 500 is weirdly high in agilent 
agilent = agilent[,c(1:499,501:512)]
par(mfrow=c(2,1))
boxplot(affy[,c(1,5,10,50,100,150,200,211,333,444,499,500,501,511)], par(las=2, cex=0.8), main='Affymetrix RMA normalised')
boxplot(agilent[,c(1,5,10,50,100,150,200,211,333,444,499,500,501,511)], par(las=2, cex=0.8), main='Agilent Lowess normalised')

# The design matrix
design = readTargets('140109_targets.txt', sep='\t')
# Remove NAs
noNAs = !is.na(design$status)
design2 = design[noNAs,]

affyDesign = overlapClinicalExpression(design2, affy)
agilentDesign = overlapClinicalExpression(design2, agilent)

# Affy was subjected gene_rma__data
affyEst = ExpressionSet(assayData=as.matrix(affyDesign[[2]]))
designAffy = affyDesign[[1]]
designAffy$status = as.factor(designAffy$status)
# May need to normalise the data at the end

# Agilent data was => unc_lowess_normalization_gene_level__data
# Fix where certsin columns have null values. This stupidly gets converted to a character matrix
aM = agilentDesign[[2]]
aM2 = data.matrix(aM)
agilentEst = ExpressionSet(assayData=aM2)
designAgilent = agilentDesign[[1]]

# Some differential expression testing. THIS currently doesn't work pick up tomorrow

# Agilent has better normalisation characteristics
dAgilent = model.matrix(~age + status, designAgilent)

# PROBLEM > designAgilent1 = !is.na(designAgilent$status)
# PROBLEM Error in design[obs, , drop = FALSE] : (subscript) logical subscript too long