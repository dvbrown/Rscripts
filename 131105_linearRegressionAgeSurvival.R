library(car)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/')
dm = read.csv('designMatrix.csv')

sur = lm(scale(dm$survival) ~ scale(dm$age))
summary(sur)
par(mfrow=c(2,2))
plot(sur)

# Assessing Outliers
par(mfrow=c(2,1))
outlierTest(sur) # Bonferonni p-value for most extreme obs
qqPlot(sur, main="QQ Plot") #qq plot for studentized resid
leveragePlots(sur) # leverage plots 

sur = lm((dm$survival) ~ (dm$age))
summary(sur)
par(mfrow=c(1,1))
plot(dm$survival ~ dm$age, main='Influence of age on survival', ylab='Survival (months)', xlab='Age (years)')
abline(sur, col ='blue')

legend("right", bty="n", legend=paste("y = -0.281x \nR2 is", format(summary(sur)$adj.r.squared, digits=3)))