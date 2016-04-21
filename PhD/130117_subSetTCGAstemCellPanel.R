library(survival)
library(limma)
source('~/Documents/Rscripts/120704-sortDataFrame.R')
source('~/Documents/Rscripts/121112_subSetTCGAfunctions.R')

#the tcga z Score transformed data
tcga = read.delim('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/121109_AgilentPooledzTransform.txt', row.names=1)
tcga = as.data.frame(tcga)
#stem cell gene signature
signature = read.delim('~/Documents/stemCellSig/130117_signature/130117_stemCellSig.txt', header=F, stringsAsFactors=F)
signature = signature[,1]
             
subsetTCGA = tcga[signature,] #subset the tcga data with the gene list. This is the signature     
upDown = (sort(apply(subsetTCGA, 1, median))) #get the median expression score to define up and down genes
#subset the gene score for upregulated and downregulated genes
up = upDown[upDown > 0] #upregulaed genes input for the geneSignature score function
down = upDown[upDown < 0]

#compute the geneScore
geneScore = betterGeneScore(tcga, up, down)
#attach clinical data to the signature score
data = bindSignatureToSurvival(geneScore)
write.table(data, '~/Documents/stemCellSig/130117_signature/130121_stemCellSig50thpercentile.txt', sep='\t', row.names=F)

#some graphs
subsetTCGA.1 = subsetTCGA[c(1:11,13,14),]
boxplot(t(subsetTCGA.1),las=2, cex.axis=1.5, main='Stem cell gene expression',cex.main=2.2, col=rainbow(13))
hist(geneScore, freq=F, col='royalblue', main='Stem cell score distribution in TCGA dataset', xlab='Stem cell signature score')
plot(data$survival, data$sigScore,main='GMB TCGA total data',ylab='Survival (days)', xlab='Stem cell signature score')
qqnorm(data$sigScore,distribution='norm')
qqline(data$sigScore, distribution=qnorm)
plot(dataSub$sigScore,dataSub$survival, main='GMB TCGA upper an lower quartiles',ylab='Survival (days)', xlab='Stem cell signature score')

#generate a column listing above the kth percentile
percentile = quantile(data$sigScore, probs=0.50, names=T)
data$percentile = ifelse(data$sigScore >= percentile, 'high', 'low')

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(data$survival, event=data$censorship)
sur.fit =survfit(data.surv~data$percentile)
plot(sur.fit, main='GBM TCGA data set',ylab='Survival probability',xlab='survival(days)', col=c('red','blue'),xlim=c(0,3000))
legend('topright', c('Stem cell score > 50th, n=220', 'Stem cell score < 50th, n=203'), col=c('red', 'blue'),lwd=1, cex=0.6)
text(locator(1),labels='p=0.0226 (Chisq)', cex=1)

summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~data$percentile)
test
