# This script performs the differential expression testing
library(affy)
library(limma)
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
list.files()

overlapClinicalExpression = function(designMatrix, expresionDataFrame) {
  # Return a list with the first element is the overlapped design matrix. The second element is the gene expression matrix
  row.names(designMatrix) = designMatrix[,1]
  dSample = designMatrix[,1]
  eSample = colnames(expresionDataFrame)
  overlap = intersect(eSample, dSample)
  dOutput = designMatrix[overlap,]
  sOutput = expresionDataFrame[,overlap]
  output = list(dOutput, sOutput)
  return (output)
}

countSurvivalStatus <- function (designMatrix) {
  # How many long and short survivors in the dataset?
  long = length(designMatrix$status[(designMatrix$status == 'long')])
  short = length(designMatrix$status[(designMatrix$status == 'short')])
  survival = list(short, long)
  names(survival) = c('short', 'long')
  return (survival)
}

affy = read.delim('140110_affyNoNulls.txt')

agilent = read.delim('140110_agilentNoNulls.txt')
# Patient number 500 is weirdly high in agilent 
agilent = agilent[,c(1:499,501:512)]
par(mfrow=c(2,1))
boxplot(affy[,c(1,5,10,50,100,150,200,211,333,444,499,500,501,511)], par(las=2, cex=0.8), main='Affymetrix RMA normalised')
boxplot(agilent[,c(1,5,10,50,100,150,200,211,333,444,499,500,501,511)], par(las=2, cex=0.8), main='Agilent Lowess normalised')

#####################################################  The design matrix ##################################################### 
design = readTargets('140113_targets_3years.txt', sep='\t')
########################################################################################################## 

# Remove NAs
noNAs = !is.na(design$status)
design2 = design[noNAs,]

affyDesign = overlapClinicalExpression(design2, affy)
agilentDesign = overlapClinicalExpression(design2, agilent)

# Affy was subjected gene_rma__data
affyEst = ExpressionSet(assayData=as.matrix(affyDesign[[2]]))
designAffy = affyDesign[[1]]
designAffy$status = as.factor(designAffy$status)

# Count the number of long and short term survivors
affyReplicates = countSurvivalStatus(designAffy)
affyReplicates
# May need to normalise the data at the end

# Agilent data was => unc_lowess_normalization_gene_level__data
# Fix where certsin columns have null values. This stupidly gets converted to a character matrix
aM = agilentDesign[[2]]
aM2 = data.matrix(aM)
agilentEst = ExpressionSet(assayData=aM2)
designAgilent = agilentDesign[[1]]

# Count the number of long and short term survivors
agilentReplicates = countSurvivalStatus(designAgilent)
agilentReplicates

##################################################### Agilent DE testing #########################################
# Agilent has better normalisation characteristics
dAgilent = model.matrix(~age + status, designAgilent)
# Fit the linear model
fitAgilent = lmFit(agilentEst, dAgilent)
# Estimate dispersions
fitAgilent = eBayes(fitAgilent)
# Specify 
resultAgilent = topTable(fitAgilent, number=17814 ,coef='statusshort', sort.by='B', adjust.method='BH')
sigAgilent = as.data.frame(decideTests(fitAgilent, p.value=0.1, lfc=1))
#write.table(resultAgilent, './limmaResults/140113_agilentShortvsLong_14months.txt', sep='\t', row.names=F)

par(mfrow=c(1,2))
volcanoplot(fitAffy, coef=3, highlight=25, names=fitAffy$genes, main='Affymetrix 3 years')

#####################################################################################################################

##################################################### AffyMetrix DE testing #########################################
dAffy = model.matrix(~age + status, designAffy)
# Fit the linear model
fitAffy = lmFit(affyEst, dAffy)
# Estimate dispersions
fitAffy = eBayes(fitAffy)
# Specify 
resultAffy = topTable(fitAffy, number=12042 ,coef='statusshort', sort.by='B', adjust.method='BH')
sigAffy = as.data.frame(decideTests(fitAffy, p.value=0.1, lfc=1))
#write.table(resultAffy, './limmaResults/140113_affymetrixShortvsLong_14months.txt', sep='\t', row.names=F)

volcanoplot(fitAgilent, coef=3, highlight=25, names=fitAgilent$genes, main='Agilent 3 years')
#####################################################################################################################

# Compare the 2 array types and the results they find
commonTable = intersect(resultAffy$ID, resultAgilent$ID)
row.names(resultAffy) = resultAffy$ID
row.names(resultAgilent) = resultAgilent$ID
agilentCommon = resultAgilent[commonTable,]
affyCommon = resultAffy[commonTable,]
mergedResult = merge.data.frame(affyCommon, agilentCommon, by.x='ID', by.y='ID')

par(mfrow=c(2,2))
plot(mergedResult$logFC.x, mergedResult$logFC.y, xlab='Affymetrix', ylab='Agilent', main='Log fold change comparsion \nbetween array platforms TCGA GBM')
plot(mergedResult$B.x, mergedResult$B.y, xlab='Affymetrix', ylab='Agilent', main='B value comparsion \nbetween array platforms TCGA GBM')

commonSig = intersect(row.names(sigAffy), row.names(sigAgilent))
agilentC = sigAgilent[commonSig,]
affyC = sigAffy[commonSig,]
plot(affyC$statusshort, agilentC$statusshort, xlab='Affymetrix', ylab='Agilent', main='Significance calls between\narray platforms TCGA GBM')
plot(affyC[,1], agilentC[,1], xlab='Affymetrix', ylab='Agilent', main='Intercept between\narray platforms TCGA GBM')