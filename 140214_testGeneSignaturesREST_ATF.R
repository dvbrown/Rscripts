# This script will take Agilent gene expression from the TCGA and attempt to subset it for survival
library(survival)
library(limma)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

############################# Check that the data has been z-transformed by looking at a few genes ############################ 
par(mfrow=c(2,2))
geneMatrix = t(agilent)
hist((geneMatrix[,'TP53']), main='p53 distribution', xlab='TP53')
hist((geneMatrix[,'GAPDH']), main='GAPDH distribution', xlab='GAPDH')
hist((geneMatrix[,'EGFR']), main='EGFR distribution', xlab='EGFR')
hist((geneMatrix[,'ACTB']), main='B-Actin distribution', xlab='B-Actin')
par(mfrow=c(1,1))

###########################################  Perform a z transformation myself ########################################## 
zScore = apply(t(geneMatrix), 1, scale)

colnames(zScore) = colnames(geneMatrix)
row.names(zScore) = row.names(geneMatrix)

par(mfrow=c(2,2))
hist((zScore[,'TP53']), main='p53 distribution', xlab='TP53')
hist((zScore[,'GAPDH']), main='GAPDH distribution', xlab='GAPDH')
hist((zScore[,'EGFR']), main='EGFR distribution', xlab='EGFR')
hist((zScore[,'ACTB']), main='B-Actin distribution', xlab='B-Actin')
par(mfrow=c(1,1))

zScore = t(zScore)
#write.table(zScore, '140218_zTransformedTCGA.txt', sep='\t')

############################################# Subset the dataframe with the signature of genes ##################################
# A dummy signature until I generate the proper REST and ATF-2 one
signature = read.delim('~/Documents/stemCellSig/130117_signature/130117_stemCellSig.txt', header=F, stringsAsFactors=F)
signature = signature[,1]

subsetTCGA = zScore[signature,] #subset the tcga data with the gene list. This is the signature     
upDown = (sort(apply(subsetTCGA, 1, median))) #get the median expression score to define up and down genes
#subset the gene score for upregulated and downregulated genes
up = upDown[upDown > 0] #upregulaed genes input for the geneSignature score function
down = upDown[upDown < 0]

#compute the geneScore
geneScore = betterGeneScore(zScore, up, down)
#attach clinical data to the signature score
data = bindSignatureToSurvival(geneScore)
write.table(data, '140220_stemCellSig50thpercentile.txt', sep='\t', row.names=F)

#some graphs
subsetTCGA.1 = subsetTCGA[c(1:11,13,14),]
boxplot(t(subsetTCGA.1),las=2, cex.axis=1.5, main='Stem cell gene expression',cex.main=2.2, col=rainbow(13))
hist(geneScore, freq=F, col='royalblue', main='Stem cell score distribution in TCGA dataset', xlab='Stem cell signature score')

plot(data$survival, data$sigScore,main='GMB TCGA total data',ylab='Survival (days)', xlab='Stem cell signature score')
qqnorm(data$sigScore,distribution='norm')
qqline(data$sigScore, distribution=qnorm)
plot(dataSub$sigScore,dataSub$survival, main='GMB TCGA upper an lower quartiles',ylab='Survival (days)', xlab='Stem cell signature score')
