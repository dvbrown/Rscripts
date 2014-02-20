# This script will take Agilent gene expression from the TCGA and attempt to subset it for survival
source('~/Documents/Rscripts/140220_TCGAsignatureAnalysisFunctions.R')

# I accidently deleted the read delim to Agilent data. Find it when I can be bothered
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/140214_testSignaturesAgilent/')

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

# Transform the clinical data into a form that can be analysed by survival package
surv = censorData(clinical)

# Compute the gene signature score
sigScore = computeSignatureScore(zScore, signature)

#attach censored clinical data to the signature score
data = bindSignatureToSurvival(sigScore, censored)
write.table(data, '140220_stemCellSignatureScore', sep='\t', row.names=F)

#some graphs to view the distrubution of the scores
subsetTCGA.1 = zScore[c(1:11,13,14),]
boxplot(t(subsetTCGA.1),las=2, cex.axis=1, main='Stem cell gene expression',cex.main=2.2, col=rainbow(13))
hist(sigScore, freq=F, col='royalblue', main='Stem cell score distribution in TCGA dataset', xlab='Stem cell signature score')

plot(data$survival, data$sigScore,main='GMB TCGA total data',xlab='Survival (days)', ylab='Stem cell signature score')
qqnorm(data$sigScore, main='Distrubtion of patients gene signature scores')
qqline(data$sigScore, distribution=qnorm, col='red')

############################################# Analysing the data for survival ##################################
