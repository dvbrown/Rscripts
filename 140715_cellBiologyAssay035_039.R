<<<<<<< HEAD 
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
=======
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

############################################## IO ###############################################
>>>>>>> 343ad01c8f7b2cbd48b9399071531aad2d73bf2b
setwd('~/Documents/Cell_biology/proliferation/Resazurin/140710_039_035/')
setwd('~/Documents/')
growthD3 = read.delim('140709_day3Rep.txt')
growthD7 = read.delim('140714_day7Rep.txt')

growthD3$clone = as.factor(growthD3$clone)
growthD7$clone = as.factor(growthD7$clone)
growthD3$subPop = as.factor(growthD3$subPop)
growthD7$subPop = as.factor(growthD7$subPop)

############################################## means and averages ###############################################

# Subtract background and take mean and SD
day3Growth = backgroundMeanSD(growthD3)
day7Growth = backgroundMeanSD(growthD7)

#plot replicates
par(mfrow=c(2,1))
plot(growthD7$rep1, growthD7$rep2, ylab='replicate 2', xlab='replicate1', main='Consistency day7', pch=16)
abline(lm(growthD7$rep1~growthD7$rep2), col='red')
plot(growthD3$rep1, growthD3$rep2, ylab='replicate 2', xlab='replicate1', main='Consistency day3', pch=16)
abline(lm(growthD3$rep1~growthD3$rep2), col='red')
=======
# Did not plate down the doublr positive
day3Growth = day3Growth[!day3Growth$subPop %in% 'CD44+/CD133+',]
day7Growth = day7Growth[!day7Growth$subPop %in% 'CD44+/CD133+',]

#### Plot the raw results ####
growthPlot3 = ggplot(data=day3Growth[day3Growth$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=subPop)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", "forestgreen", "red")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 3 by \nmarker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

growthPlot7 = ggplot(data=day7Growth[day7Growth$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=subPop)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", "forestgreen")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 \nby marker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#multiplot(growthPlot3, growthPlot7)

############################################## Calculate and plot the DMSO corrected values #################################################
day3TMZ = calcDMSOcontrol(day3Growth)
day7TMZ = calcDMSOcontrol(day7Growth)
# 
tmzPlot3 = ggplot(data=day3TMZ, aes(x=clone, y=dmsoCorrected, fill=subPop)) + 
    scale_fill_manual(values=c("gold", "chartreuse4", "skyblue2")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Clone") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Temozolomide sensitivty at day 3 by \nmarker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

tmzPlot7 = ggplot(data=day7TMZ, aes(x=clone, y=dmsoCorrected, fill=subPop)) + 
    scale_fill_manual(values=c("gold", "chartreuse4", "skyblue2")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Clone") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Temozolomide sensitivty at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#multiplot(tmzPlot3, tmzPlot7)

#multiplot(growthPlot3, growthPlot7, tmzPlot3, tmzPlot7, cols=2)
# write.table(day3TMZ, '140715_day3TMZprocessed.txt', sep='\t')
# write.table(day7TMZ, '140715_day7TMZprocessed.txt', sep='\t')
>>>>>>> 343ad01c8f7b2cbd48b9399071531aad2d73bf2b
