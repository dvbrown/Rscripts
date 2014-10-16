source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data and then remove background fluorescence and take the mean and sd
    # Subtract the background (0 cells) which is the first row
    background = rowMeans(dataFrame[1,4:6], na.rm=T)
    # Subtract the background and bind back the metaData
    holder = dataFrame[,c(4:6)] - background
    result = cbind(dataFrame[,c(7:9)], holder)
    # Compute the mean and sd
    result$mean = rowMeans(result[,c(4:6)], na.rm=T)
    result$sd = apply(result[,c(4:6)], 1, sd, na.rm=T)
    return (result)
}

calcDMSOcontrol = function(dataFrame) {
    vehicle = dataFrame[dataFrame$treatment %in% 'DMSO',]
    tmz = dataFrame[dataFrame$treatment %in% 'TMZ',]
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

summariseByFactor <- function (dataFrame, factor1, factor2) {
    # Takes a dataframe with factor information and computes summary statistics based on levels of as least 2 factors
    # factor 1 and 2 are characters
    result <- ddply(dataFrame, c(factor1, factor2), summarise,
                   N    = length(percent),
                   mean = mean(percent),
                   sd   = sd(percent),
                   se   = sd / sqrt(N) )
  return (result)
}

takeUnison <- function (clinicalData, genomicData) {
    # Takes the union of cases in clinical and genomic data.
    # Sort both datasets so the order is the same
    samples = intersect(colnames(genomicData), row.names(clinicalData))
    clinIntersect = clinicalData[samples,]
    gemIntersect = genomicData[,samples]
    clinIntersect <- clinIntersect[order(row.names(clinIntersect)),]
    gemIntersect <- gemIntersect[order(row.names(gemIntersect)),]
    # Returns a list with clinical data first and genomic data second
    result = list(clinIntersect, gemIntersect)
}