#row bind all expression data and z-score transform
setwd('~/Documents/public-datasets/TCGA/expressionData/')

zTransform = function(matrixElement, rowMean, rowSD ) { #convert to the unitless z-score based on a normal distribution
    z = (matrixElement - rowMean)/rowSD
    return (z)
}

agilent1 = read.delim('Expression-Genes/UNC__AgilentG4502A_07_1/120727-Aglient1SummarisedData.txt', row.names=1)
agilent2 = read.delim('Expression-Genes/UNC__AgilentG4502A_07_2/121106_Agilent2data.txt', row.names=1)
affy = read.delim('Expression-Genes/BI__HT_HG-U133A/121107_collatedExpressionAffy.txt', row.names=1)

agilentBind = cbind(agilent1, agilent2)
rowMean = rowMeans(agilentBind, na.rm=T)
dataMatrix = apply(agilentBind, 2, as.numeric)
row.names(dataMatrix) = row.names(agilentBind)
rowMean = rowMeans(dataMatrix)
rowStdDev = apply(dataMatrix, 1, sd)

zScore = apply(dataMatrix, 2, zTransform, rowMean, rowStdDev) #compute the z-scores for the dataFrame