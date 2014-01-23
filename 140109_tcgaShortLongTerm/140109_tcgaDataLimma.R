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

#####################################################  The design matrix #################################
# Female is '1'. Male is '2'
design = readTargets('140115_design3yearMatrix.txt', sep='\t')
########################################################################################################## 

affy = read.delim('140110_affyNoNulls.txt')

agilent = read.delim('140110_agilentNoNulls.txt')
# Patient number 500 is weirdly high in agilent 
agilent = agilent[,c(1:499,501:512)]
par(mfrow=c(2,1))
boxplot(affy[,c(1,5,10,50,100,150,200,211,333,444,499,500,501,511)], par(las=2, cex=0.8), main='Affymetrix RMA normalised')
boxplot(agilent[,c(1,5,10,50,100,150,200,211,333,444,499,500,501,511)], par(las=2, cex=0.8), main='Agilent Lowess normalised')

# Remove NAs
noNAs = !is.na(design$status)
design2 = design[noNAs,]

affyDesign = overlapClinicalExpression(design2, affy)
agilentDesign = overlapClinicalExpression(design2, agilent)

# Affy was subjected gene_rma__data
aF = affyDesign[[2]]
aF2 = data.matrix(aF)

########################################### Normalise Affymetrix #############################################
# It is necesary to log transform the Affymetrix data
ks.test(aF2[400,], 'pnorm', mean=mean(aF2[400,]), sd=sd(aF2[400,]))
aF3 = log2(aF2)
aF3 = normalizeBetweenArrays(aF3, method='quantile')
# Some plots to explore distribution
ks.test(aF3[400,], 'pnorm', mean=mean(aF3[400,]), sd=sd(aF3[400,]))
boxplot(aF3[,c(1,5,10,50,100,150,200,211,333,375,400,409)], main='Affymetrix quantile normalised')#, par(las=2, cex=0.8))
hist(aF2[,1], breaks='FD', main='Affymetrix patient 1 expression', xlab='RMA normalised expression score')
qqnorm(aF2[,1], main='Affymetrix patient 1 Q-Q')
qqline(aF2[,1])

gapdhAf = aF2['ACTB',]
hist(gapdhAf, breaks='FD', main='Affymetrix B-Actin \nexpression distribution')
qqnorm(gapdhAf, main='Affymetrix B-Actin Q-Q')
qqline(gapdhAf)

affyEst = ExpressionSet(assayData=aF3)
designAffy = affyDesign[[1]]
designAffy$status = as.factor(designAffy$status)

# Count the number of long and short term survivors
affyReplicates = countSurvivalStatus(designAffy)
affyReplicates
##############################################################################################################

# Agilent data was => unc_lowess_normalization_gene_level__data
aM = agilentDesign[[2]]
aM2 = data.matrix(aM)
#boxplot(aM3[,c(1,5,10,50,100,150,200,211,333,395)], main='Agilent Lowess + Scale normalised')#, par(las=2, cex=0.8))

# Some plots to explore distribution
hist(aM2[,1], breaks='FD', main='Agilent patient 1 expression', xlab='Loess normalised expression score')
ks.test(aM2[250,], 'pnorm')#, mean=mean(aM3[,300]), sd=sd(aM3[,300]))
qqnorm(aM2[,c(1)],main='Agilent patient 1 Q-Q')
qqline(aM2[,1])

gapdhAg = aM2['ACTB',]
hist(gapdhAg, breaks='FD', main='Agilent B-Actin \nexpression distribution')
qqnorm(gapdhAg, main='Agilent B-Actin Q-Q')
qqline(gapdhAg)

agilentEst = ExpressionSet(assayData=aM2)
designAgilent = agilentDesign[[1]]

# Count the number of long and short term survivors
agilentReplicates = countSurvivalStatus(designAgilent)
agilentReplicates

##################################################### Agilent DE testing #########################################
# Agilent has better normalisation characteristics
dAgilent = model.matrix(~age + gender + status, designAgilent)

# Fit the linear model
fitAgilent = lmFit(agilentEst, dAgilent)
# Estimate dispersions
fitAgilent = eBayes(fitAgilent)
# Specify 
resultAgilent = topTable(fitAgilent, number=17814 ,coef='statusshort', sort.by='B', adjust.method='BH')
sigAgilent = as.data.frame(decideTests(fitAgilent, p.value=0.1, lfc=1))
#write.table(resultAgilent, './limmaResults/140115_agilentShortvsLong_3yearsGender.txt', sep='\t', row.names=F)

##################################################### AffyMetrix DE testing #########################################
dAffy = model.matrix(~age + gender + status, designAffy)
# Fit the linear model
fitAffy = lmFit(affyEst, dAffy)
# Estimate dispersions
fitAffy = eBayes(fitAffy)
# Specify 
resultAffy = topTable(fitAffy, number=12042 ,coef='statusshort', sort.by='B', adjust.method='BH')
sigAffy = as.data.frame(decideTests(fitAffy, p.value=0.1, lfc=1))
#write.table(resultAffy, './limmaResults/140115_affymetrixShortvsLong_3yearsGender.txt', sep='\t', row.names=F)

# Draw the volcano plots for both Agilent and Affymetrix
#par(mfrow=c(1,2))
#volcanoplot(fitAffy, coef=3, highlight=25, names=fitAffy$genes, main='Affymetrix 3 years')
#legend('topleft', 'short-term n=378 \nlong-term=31', cex=0.75)
#volcanoplot(fitAgilent, coef=3, highlight=25, names=fitAgilent$genes, main='Agilent 3 years')
#legend('topleft', 'short-term n=363 \nlong-term=32', cex=0.75)
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