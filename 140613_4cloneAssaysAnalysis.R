source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data and then remove background fluorescence and take the mean and sd
    # Subtract the background (0 cells) which is the first row
    background = rowMeans(dataFrame[21,4:6], na.rm=T)
    # Subtract the background and bind back the metaData
    holder = dataFrame[,c(4:6)] - background
    result = cbind(dataFrame[,c(7:9)], holder)
    # Compute the mean and sd
    result$mean = rowMeans(result[,c(4:6)], na.rm=T)
    result$sd = apply(result[,c(4:6)], 1, sd, na.rm=T)
    result$treatment = dataFrame$treatment
    return (result)
}

calcDMSOcontrol = function(dataFrame) {
    vehicle = dataFrame[dataFrame$treatment %in% 'DMSO',]
    tmz = dataFrame[dataFrame$treatment %in% 'TMZ',]
    tmz$dmsoCorrected = tmz$mean / vehicle$mean
    return (tmz)
}

calcProlifNormalised = function(dataFrame) {
    negative = dataFrame[dataFrame$cd133 %in% 'CD133_neg',]
    positive = dataFrame[dataFrame$cd133 %in% 'CD133_pos',]
    positive$negNormalised = positive$mean / negative$mean
    # Now the double stains
    doubleNeg = dataFrame[dataFrame$cd133 %in% 'doubleNeg',]
    otherPop = dataFrame[!dataFrame$cd133 %in% c('doubleNeg','CD133_neg','CD133_pos'),] 
    otherPop$negNormalised = otherPop$mean / doubleNeg$mean
    result = rbind(positive, otherPop)
    return (result)
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
setwd('~/Documents/Cell_biology/proliferation/Resazurin/140602_4clones/inputs/analysis/')
day7 = read.delim('../140613_day7RinputV2.txt')

plot(day7$rep1, day7$rep2, ylab='replicate 2', xlab='replicate1', main='Consistency day7', pch=16, col=day7$clone)
abline(lm(day7$rep1~day7$rep2), col='red')

day7Growth = backgroundMeanSD(day7)
day7Growth = day7Growth[c(1:20),]
#write.table(day7Growth, '140613_day7meanSD.txt', sep='\t')

# Plot the raw results of growth
growthPlot7 = ggplot(data=day7Growth[day7Growth$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=cd133status)) + 
    #scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing proliferation at day 7 \nby CD133/ CD44 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

tmzRaw7 = ggplot(data=day7Growth[day7Growth$treatment %in% 'TMZ',], 
                 aes(x=clone, y=mean, fill=cd133status)) + 
    #scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing TMZ at day 7 by \nCD133/ CD44 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(growthPlot7, tmzRaw7)

############################################## Calculate the DMSO corrected values #################################################
day7TMZ = calcDMSOcontrol(day7Growth)

tmzPlot7 = ggplot(data=day7TMZ, aes(x=clone, y=dmsoCorrected, fill=cd133status)) + 
    #scale_fill_manual(values=c("gold", "chartreuse4")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Clone") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Comparing temozolomide sensitivty at day 7\nby CD133/CD44 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#write.table(day7TMZ, '140613_day7TMZprocessed.txt', sep='\t')

############################################## Normalise proliferation to CD133 negative #################################################
day7GrowthNorm = calcProlifNormalised(day7Growth[day7Growth$treatment %in% 'DMSO',])
day7GrowthNorm$ID = paste(day7GrowthNorm$clone, day7GrowthNorm$cd133status)

day7GrowthNormP = ggplot(data=day7GrowthNorm, aes(x=ID, y=negNormalised, fill=clone)) + 
    scale_fill_manual(values=c("darkorange", "royalblue","cyan","magenta")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Clone") + ylab("Cell number relative to matched CD133 negative") +
    ggtitle("Proliferation of CD133/ CD44 cells at day 7") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

############################################## Summarise data by FACS status #################################################
