source('140211_multiplotGgplot2.R')

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data and then remove background fluorescence and take the mean and sd
    # Subtract the background (0 cells) which is the first row
    background = rowMeans(dataFrame[4,4:6], na.rm=T)
    # Subtract the background and bind back the metaData
    holder = dataFrame[,c(4:6)] - background
    result = cbind(dataFrame[,c(7:9)], holder)
    # Compute the mean and sd
    result$mean = rowMeans(result[,c(4:6)], na.rm=T)
    result$sd = apply(result[,c(4:6)], 1, sd, na.rm=T)
    return (result)
}

calcDMSOcontrol = function(dataFrame) {
    vehicle = dataFrame[dataFrame$treatment %in% 'vehicle',]
    tmz = dataFrame[dataFrame$treatment %in% 'tmz',]
    tmz$dmsoCorrected = tmz$mean / vehicle$mean
    return (tmz)
}

calcProlifNormalised = function(dataFrame) {
    negative = dataFrame[dataFrame$cd133 %in% 'neg',]
    positive = dataFrame[dataFrame$cd133 %in% 'pos',]
    positive$negNormalised = positive$mean / negative$mean
    return (positive)
}

extractPosNegReplicates = function(dataFrame) {
    neg = dataFrame[dataFrame$cd133 %in% 'neg',c(1,2,3,9)]
    pos = dataFrame[dataFrame$cd133 %in% 'pos',c(1,2,3,9)]
    negMean = mean(neg$dmsoCorrected)
    negSD = sd(neg$dmsoCorrected)
    posMean = mean(pos$dmsoCorrected)
    posSD = sd(pos$dmsoCorrected)
    negSummary = c(negMean, negSD)
    posSummary = c(posMean, posSD)
    result = rbind(negSummary, posSummary)
    result = as.data.frame(result)
    origin = c('negative', 'positive')
    result = cbind(origin, result)
    colnames(result) = c('origin', 'mean', 'sd')
    return (result)
}

############################################## Read in the resazurin assay readings ###############################################
setwd('~/Documents/Cell_biology/proliferation/Resazurin/140710_039_035/')
setwd('~/Documents/')
growthD3 = read.delim('140709_day3Rep.txt')
growthD7 = read.delim('140714_day7Rep.txt')

# Subtract background and take mean and SD
day3Growth = backgroundMeanSD(growthD3)
day7Growth = backgroundMeanSD(growthD7)

#plot replicates
par(mfrow=c(2,1))
plot(growthD7$rep1, growthD7$rep2, ylab='replicate 2', xlab='replicate1', main='Consistency day7', pch=16)
abline(lm(growthD7$rep1~growthD7$rep2), col='red')
plot(growthD3$rep1, growthD3$rep2, ylab='replicate 2', xlab='replicate1', main='Consistency day3', pch=16)
abline(lm(growthD3$rep1~growthD3$rep2), col='red')