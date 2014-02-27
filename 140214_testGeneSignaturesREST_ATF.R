# This script will take Agilent gene expression from the TCGA and attempt to subset it for survival
source('~/Documents/Rscripts/140220_TCGAsignatureAnalysisFunctions.R')

agilent = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/140110_affyNoNulls.txt')
# I accidently deleted the read delim to Agilent data. Find it when I can be bothered
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/140214_testSignaturesAgilent/')

############################# Check that the data has been z-transformed by looking at a few genes ############################ 
# par(mfrow=c(2,2))
# geneMatrix = t(agilent)
# hist((geneMatrix[,'TP53']), main='p53 distribution', xlab='TP53')
# hist((geneMatrix[,'GAPDH']), main='GAPDH distribution', xlab='GAPDH')
# hist((geneMatrix[,'EGFR']), main='EGFR distribution', xlab='EGFR')
# hist((geneMatrix[,'ACTB']), main='B-Actin distribution', xlab='B-Actin')
# par(mfrow=c(1,1))
# 
# ###########################################  Perform a z transformation myself ########################################## 
# zScore = apply(t(geneMatrix), 1, scale)
# 
# colnames(zScore) = colnames(geneMatrix)
# row.names(zScore) = row.names(geneMatrix)
# 
# par(mfrow=c(2,2))
# hist((zScore[,'TP53']), main='p53 distribution', xlab='TP53')
# hist((zScore[,'GAPDH']), main='GAPDH distribution', xlab='GAPDH')
# hist((zScore[,'EGFR']), main='EGFR distribution', xlab='EGFR')
# hist((zScore[,'ACTB']), main='B-Actin distribution', xlab='B-Actin')
# par(mfrow=c(1,1))
# 
# zScore = t(zScore)
# #write.table(zScore, '140218_zTransformedTCGA.txt', sep='\t')

############################################# Subset the dataframe with the signature of genes ##################################
# A dummy signature until I generate the proper REST and ATF-2 one
signature = read.delim('~/Documents/RNAdata/danBatch1/GSEA/signaturesToTest/PID_ATF2_PATHWAY.txt', header=F, stringsAsFactors=F)
signature = signature[,1]

# Transform the clinical data into a form that can be analysed by survival package
surv = censorData(clinical)

# Compute the gene signature score
sigScore = computeSignatureScore(zScore, signature)

#attach censored clinical data to the signature score
data = bindSignatureToSurvival(sigScore, surv)

write.table(data, 'output.txt', sep='\t', row.names=T)
# For whatever reason censorship column is not nicley coerced to integer
rm(data)
data = read.delim('output.txt', row.names=1)
# data = data[complete.cases(data$survival),]
# data = data[complete.cases(data$censorship),]
data$censorship = as.integer(data$censorship)

#some graphs to view the distrubution of the scores
subsetTCGA.1 = zScore[c(1:11,13,14),]
boxplot(t(subsetTCGA.1),las=2, cex.axis=1, main='Stem cell gene expression',cex.main=2.2, col=rainbow(13))
hist(sigScore, freq=F, col='royalblue', main='Stem cell score distribution in TCGA dataset', xlab='Stem cell signature score')
plot(data$survival, data$sigScore,main='GMB TCGA total data',xlab='Survival (days)', ylab='Stem cell signature score')
qqnorm(data$sigScore, main='Distrubtion of patients gene signature scores')
qqline(data$sigScore, distribution=qnorm, col='red')

############################################# Analysing the data for survival ##################################

#generate a column listing above the kth percentile
data = buildClassifier(data, 0.5)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(data$survival, event=data$censorship)
sur.fit = survfit(data.surv~data$percentile)

plot(sur.fit, main='Predictive power of ATF2 signature \nin the TCGA',ylab='Survival probability',xlab='survival(days)', col=c('red','blue'),xlim=c(0,750))
legend('topright', c('Gene signature score > 50th, n=188', 'Gene signature score < 50th, n=193'), col=c('red', 'blue'),lwd=1, cex=0.6)
summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~data$percentile)
test
text(locator(1),labels='p=0.536', cex=1) #add the p-value to the graph
