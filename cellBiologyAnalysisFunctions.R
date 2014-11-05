source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

backgroundMeanSD <- function (dataFrame, blankEntry) {
    # Take the dataframe of raw data and then remove background fluorescence and take the mean and sd
    # Subtract the background which is indexed by the row number in the argument to the function
    background = rowMeans(dataFrame[blankEntry,4:6], na.rm=T)
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

# calcGrowthNormCD133 = function(dataFrame) {
#     negative = dataFrame[dataFrame$cd133 %in% 'neg',]
#     positive = dataFrame[dataFrame$cd133 %in% 'pos',]
#     positive$negNormalised = positive$mean / negative$mean
#     return (positive)
# }

calcGrowthNormalised = function(dataFrame, patientName) {
    # Normalises a patient by the double negative subpopulation set to 1
    # Patient is a charaacter string of patient eg #035. Treat is either DMSO or TMZ
    # Extract all cases of the individual patient
    patient = dataFrame[dataFrame[,"patient"] %in% patientName,]
    patient = patient[patient[,"treatment"] %in% "DMSO",]
    # Extract the double negative
    dn = patient[patient[,"subpop"] %in% "CD44-/CD133-",]
    otherSample = patient
    #     otherSample$norm1 = otherSample$rep1 / dn$rep1
    #     otherSample$norm3 = otherSample$rep3 / dn$rep3
    #     otherSample$norm2 = otherSample$rep2 / dn$rep2
    otherSample$normDN = otherSample$mean / dn$mean
    return (otherSample)
}

calcTMZNormalised = function(dataFrame, patientName) {
    # Normalises a patient by the double negative subpopulation set to 1
    # Patient is a charaacter string of patient eg #035. Treat is either DMSO or TMZ
    # Extract all cases of the individual patient
    patient = dataFrame[dataFrame[,"patient"] %in% patientName,]
    # Extract the double negative
    dn = patient[patient[,"subpop"] %in% "CD44-/CD133-",]
    otherSample = patient
    otherSample$normDN = otherSample$dmsoCorrected / dn$dmsoCorrected
    return (otherSample)
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
    require(plyr)
    # Takes a dataframe with factor information and computes summary statistics based on levels of as least 2 factors
    # factor 1 and 2 are characters
    result <- ddply(dataFrame, c(factor1, factor2), summarise,
                   N    = length(normDN),
                   mean = mean(normDN),
                   sd   = sd(normDN),
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