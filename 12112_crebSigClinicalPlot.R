library(survival)
library(limma)

crebSigScore = read.delim('./Documents/CREB/CREBsignatures/121112_MontminyPromoter_zScore/121112_paitentSignature.txt',row.names=1)
clinicalData = read.delim('~/Documents/public-datasets/TCGA/clinicalData/120731-survivalDataStructered.txt', row.names=1)
clinSub = clinicalData[,c(4,5)]
finalData = merge.data.frame(clinSub, crebSigScore, by.x='row.names', by.y='row.names')
plot(finalData$censor, finalData$crebScore, main='CRE containing promoters and survival', ylab='CREB signature score',
     xlab='survival')
#Create a new vector based on CREB score. If below 0 code in 1, otherwise code in 2. Change this for grouping
finalData$group = ifelse(finalData$crebScore < 0, 1,2 )

subsetQuantiles = function(sigSurvivalDf, lowerQuartile, upperQuartile) { #input is a integer for the quantiles
  r = quantile(finalData$crebScore, c(lowerQuartile, upperQuartile)) #get the quantiles for creb signature score. Use this to set group above
  cutFinalDataHi = subset.data.frame(finalData, crebScore > r[2])
  cutFinalDataLo = subset.data.frame(finalData, crebScore < r[1])
  data = rbind(cutFinalDataLo, cutFinalDataHi)
  return(data)
}

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(finalData$censor, event=finalData$status)
sur.fit =survfit(data.surv~group, data=finalData) #code in the quartiles
plot(sur.fit, main='Kaplan-Meier high and low CREB values',ylab='probability',xlab='survival(days)', col=c('red','blue'),xlim=c(0,3000))
legend('topright', c('25th percentile CREB score', '75th percentile CREB score'), col=c('blue', 'red'),lwd=1)
summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~surData$gender)