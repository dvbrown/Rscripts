setwd('~/Documents/Cell_biology/proliferation/Resazurin/140614_summary/analysisByMarkerBias/')
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

calcDMSOcontrol = function(dataFrame) {
    vehicle = dataFrame[dataFrame$treatment %in% 'DMSO',]
    tmz = dataFrame[dataFrame$treatment %in% 'TMZ',]  
    tmz$rep1 = tmz$rep1 / vehicle$rep1
    tmz$rep2 = tmz$rep2 / vehicle$rep2
    tmz$rep3 = tmz$rep3 / vehicle$rep3
    tmz$mean = rowMeans(tmz[,c(4:6)], na.rm=T)
    tmz$sd = apply(tmz[,c(4:6)], 1, sd, na.rm=T)
    return (tmz)
}

calcProlifNormalised = function(dataFrame) {
    negative = dataFrame[dataFrame$cd133 %in% 'CD133_neg',]
    positive = dataFrame[dataFrame$cd133 %in% 'CD133_pos',]
    positive$rep1 = positive$rep1 / negative$rep1
    positive$rep2 = positive$rep2 / negative$rep2
    positive$rep3 = positive$rep3 / negative$rep3
    positive$mean = rowMeans(positive[,c(4:6)], na.rm=T)
    positive$sd = apply(positive[,c(4:6)], 1, sd, na.rm=T)
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

normaliseCD133 <- function (dataFrame) {
    cd133Neg = dataFrame[dataFrame$cd133status %in% 'CD133_neg',]
    cd133Pos = dataFrame[dataFrame$cd133status %in% 'CD133_pos',]
    cd133NegAv = mean(cd133Neg$mean)
    cd133NegSd = sd(cd133Neg$mean) / sqrt(length(cd133Neg$mean))
    cd133PosAv = mean(cd133Pos$mean)
    cd133PosSd = sd(cd133Pos$mean) / sqrt(length(cd133Pos$mean))
    cd133 = as.data.frame(rbind(c(cd133NegAv, cd133NegSd), c(cd133PosAv, cd133PosSd)))
    cd133$cd133 = c('negative', 'positive')
    return (cd133)
}

################################ IO and subsetting ############################################
resazurin = read.delim('~/Documents/Cell_biology/proliferation/Resazurin/140614_summary/140614_day7meanSD.txt')
invasion = read.delim('~/Documents/Cell_biology/microscopy/invasion/140615_summary/140617_summary.rep.txt')
clinical = read.delim('~/Documents/Cell_biology/140617_cloneProgressAssays.txt')
invasion$clone = c('004','004','004','004','020','020','020','034','034','034',
                '034','035','035','035','035','039','039','039','039','041','041','041', '039')
invasion$clone = as.factor(invasion$clone)

# Add the population bias to the data
resazurinC = merge.data.frame(resazurin, clinical, by.x='clone', by.y='clone')
invasionC = merge.data.frame(invasion, clinical, by.x='clone', by.y='clone')
