source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

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
setwd('~/Documents/Cell_biology/proliferation/Resazurin/140417_6clones/analysis/')
growthD3 = read.delim('140414_day3_linearRep.txt')
growthD7 = read.delim('140417_day7_linearRep.txt')

# par(mfrow=c(2,1))
# plot(growthD7$rep1, growthD7$rep2, ylab='replicate 2', xlab='replicate1', main='Consistency day7', pch=16)
# abline(lm(growthD7$rep1~growthD7$rep2), col='red')
# plot(growthD3$rep1, growthD3$rep2, ylab='replicate 2', xlab='replicate1', main='Consistency day3', pch=16)
# abline(lm(growthD3$rep1~growthD3$rep2), col='red')

# Subtract background and take mean and SD
day3Growth = backgroundMeanSD(growthD3)
day7Growth = backgroundMeanSD(growthD7)

# Plot the raw results
# growthPlot3 = ggplot(data=day3Growth[day3Growth$treatment %in% 'growth',], 
#                    aes(x=clone, y=mean, fill=cd133)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
#     xlab("Clone") + ylab("Fluorescent intensity") +
#     ggtitle("Comparing proliferation at day 3 by CD133 status") +  # Set title
#     theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# growthPlot7 = ggplot(data=day7Growth[day7Growth$treatment %in% 'growth',], 
#                      aes(x=clone, y=mean, fill=cd133)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
#     xlab("Clone") + ylab("Fluorescent intensity") +
#     ggtitle("Comparing proliferation at day 7 by CD133 status") +  # Set title
#     theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# multiplot(growthPlot3, growthPlot7)

####################################################################################################################################

# tmzRaw3 = ggplot(data=day3Growth[day3Growth$treatment %in% 'tmz',], 
#                      aes(x=clone, y=mean, fill=cd133)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
#     xlab("Clone") + ylab("Fluorescent intensity") +
#     ggtitle("Comparing TMZ at day 3 by CD133 status") +  # Set title
#     theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# tmzRaw7 = ggplot(data=day7Growth[day7Growth$treatment %in% 'tmz',], 
#                  aes(x=clone, y=mean, fill=cd133)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
#     xlab("Clone") + ylab("Fluorescent intensity") +
#     ggtitle("Comparing TMZ at day 7 by CD133 status") +  # Set title
#     theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# multiplot(tmzRaw3, tmzRaw7)

############################################## Calculate the DMSO corrected values #################################################
day3GrowthMatched = day3Growth[!day3Growth$clone %in% c('030a_pos', '034a_neg', 'blank', 'empty'),]
day7GrowthMatched = day7Growth[!day7Growth$clone %in% c('030a_pos', '034a_neg', 'blank', 'empty'),]
day3TMZ = calcDMSOcontrol(day3GrowthMatched)
day7TMZ = calcDMSOcontrol(day7GrowthMatched)

# tmzPlot3 = ggplot(data=day3TMZ, aes(x=clone, y=dmsoCorrected, fill=cd133)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     xlab("Clone") + ylab("Cell number relative to DMSO control") +
#     ggtitle("Comparing temozolomide sensitivty at day 3 by CD133 status") +  # Set title
#     theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# tmzPlot7 = ggplot(data=day7TMZ, aes(x=clone, y=dmsoCorrected, fill=cd133)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     xlab("Clone") + ylab("Cell number relative to DMSO control") +
#     ggtitle("Comparing temozolomide sensitivty at day 7 by CD133 status") +  # Set title
#     theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# multiplot(tmzPlot3, tmzPlot7)

# multiplot(growthPlot3, growthPlot7, tmzPlot3, tmzPlot7)
# write.table(day3TMZ, '140423_day3TMZprocessed.txt', sep='\t')
# write.table(day7TMZ, '140423_day7TMZprocessed.txt', sep='\t')

############################################## Normalise proliferation to CD133 negative #################################################

day3GrowthNorm = calcProlifNormalised(day3GrowthMatched[day3GrowthMatched$treatment %in% 'growth',])
day7GrowthNorm = calcProlifNormalised(day7GrowthMatched[day7GrowthMatched$treatment %in% 'growth',])

# day3GrowthNormP = ggplot(data=day3GrowthNorm, aes(x=clone, y=negNormalised, fill=clone)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue","cyan","magenta")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     xlab("Clone") + ylab("Cell number relative to matched CD133 negative") +
#     ggtitle("Comparing proliferation of CD133 cells at day 3") +  # Set title
#     theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# day7GrowthNormP = ggplot(data=day7GrowthNorm, aes(x=clone, y=negNormalised, fill=clone)) + 
#     scale_fill_manual(values=c("darkorange", "royalblue","cyan","magenta")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     xlab("Clone") + ylab("Cell number relative to matched CD133 negative") +
#     ggtitle("Comparing proliferation of CD133 cells at day 7") +  # Set title
#     theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# multiplot(day3GrowthNormP, day7GrowthNormP)

# Summarise data 
dayG3sum = c(mean(day3GrowthNorm$negNormalised), sd(day3GrowthNorm$negNormalised))
dayG3sum = rbind(dayG3sum, c(1, 0))
dayG3sum = as.data.frame(dayG3sum)

dayG7sum = c(mean(day7GrowthNorm$negNormalised), sd(day7GrowthNorm$negNormalised))
dayG7sum = rbind(dayG7sum, c(1, 0))
dayG7sum = as.data.frame(dayG7sum)

growthSummary = rbind(dayG3sum, dayG7sum)
colnames(growthSummary) = c('average', 'stdDev')
row.names(growthSummary) = c(1,2,3,4)
growthSummary$origin = c('pos', 'neg', 'pos', 'neg')
growthSummary$day = c('day 3', 'day 3', 'day 7', 'day 7')

# growthSummaryP = ggplot(data=growthSummary, aes(x=day, y=average, fill=origin)) + 
#     scale_fill_manual(values=c("aquamarine1", "darkgoldenrod1")) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     geom_errorbar(aes(ymin=average-stdDev, ymax=average+stdDev), width=.2, position=position_dodge(0.9)) +
#     xlab("measurment day") + ylab("Normalised cell number") +
#     ggtitle("Growth rate of GIC clones (n = 4)") +  # Set title
#     theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# growthSummaryP


############################################## Pool the pos and negs for stats #################################################
tmzSummary3 = extractPosNegReplicates(day3TMZ)
tmzSummary7 = extractPosNegReplicates(day7TMZ)
# Sumamrise data
tmzSummary = rbind(tmzSummary3, tmzSummary7)
tmzSummary$day = c('day 3', 'day 3', 'day 7', 'day 7')

# tmzSummary3P = ggplot(data=tmzSummary3, aes(x=origin, y=mean, fill=origin)) + 
#     scale_fill_manual(values=c("aquamarine1", "darkgoldenrod1")) +
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     xlab("") + ylab("Cell number relative to DMSO control") +
#     ggtitle("Effect of temozolomide on subpopulations of GICs at day 3") +  # Set title
#     theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# tmzSummary7P = ggplot(data=tmzSummary7, aes(x=origin, y=mean, fill=origin)) + 
#     scale_fill_manual(values=c("darkcyan", "coral3")) +
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     xlab("") + ylab("Cell number relative to DMSO control") +
#     ggtitle("Effect of temozolomide on subpopulations of GICs at day 7") +  # Set title
#     theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# tmzSummary7P + annotate("text", x = 25, y = 25, label = "n = 4")
# 
# multiplot(tmzSummary3P, tmzSummary7P)

# tmzSummaryP = ggplot(data=tmzSummary, aes(x=day, y=mean, fill=origin)) + 
#     scale_fill_manual(values=c("aquamarine1", "darkgoldenrod1")) +
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
#     geom_bar(stat="identity", position=position_dodge(), colour="black") + 
#     xlab("") + ylab("Cell number relative to DMSO control") +
#     ggtitle("Effect of temozolomide on subpopulations of GICs") +  # Set title
#     theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# tmzSummaryP

####################################################################################################################################





############################################## Read in the invasion assay readings #################################################
# Clear memory
rm(list=ls())
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
setwd('~/Documents/Cell_biology/microscopy/invasion/140414_invasion/')
invD3 = read.delim('140422_outputDay3.rep.txt')
invD7 = read.delim('140422_outputDay7.rep.txt')

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data take the mean and sd
    dataFrame$mean = rowMeans(dataFrame[,c(4:6)], na.rm=T)
    dataFrame$sd = apply(dataFrame[,c(4:6)], 1, sd, na.rm=T)
    return (dataFrame)
}

normaliseMatrixCD133 = function(dataFrame) {
    # Normalise Matrix first
    noMatrix = dataFrame[dataFrame$matrix %in%  FALSE,]
    matrix = dataFrame[dataFrame$matrix %in% TRUE,]
    matrix$matNormalised = matrix$mean / noMatrix$mean
    matrix$matNormalisedSD = matrix$sd / noMatrix$sd
    # Normalise CD133
    negative = matrix[matrix$cd133 %in% 'neg',]
    positive = matrix[matrix$cd133 %in% 'pos',]
    positive$cd133Norm = positive$mean / negative$mean
    # Return both dataframes
    result = list(matrix, positive)
    return (result)
}

invD3 = invD3[!invD3$matrix %in% NA,]
invD3$sample = paste(invD3$clone, invD3$cd133, sep='_')
invD7 = invD7[!invD7$matrix %in% NA,]
invD7$sample = paste(invD7$clone, invD7$cd133, sep='_')

invD3Stats = backgroundMeanSD(invD3)
invD7Stats = backgroundMeanSD(invD7)

invD3StatsP = ggplot(data=invD3Stats[!invD3Stats$clone %in% c('030a', '034a'),], aes(x=sample, y=mean, fill=matrix)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Invasive ability of GIC clones buy CD133 status at day 3") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

invD7StatsP = ggplot(data=invD7Stats[!invD7Stats$clone %in% c('030a', '034a'),], aes(x=sample, y=mean, fill=matrix)) + 
    scale_fill_manual(values=c("yellow", "skyblue3")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Invasive ability of GIC clones buy CD133 status at day 7") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(invD3StatsP, invD7StatsP)
####################################################################################################################################

############################################## Normalise for no matrix and CD133 #################################################

# Remove the unmatched recurrents
invD3Stats = invD3Stats[!invD3Stats$clone %in% c('030a', '034a'),]
invD7Stats = invD7Stats[!invD7Stats$clone %in% c('030a', '034a'),]

invD3Norm = normaliseMatrixCD133(invD3Stats)
invD7Norm = normaliseMatrixCD133(invD7Stats)

invD3NormP = ggplot(data=invD3Norm[[1]], aes(x=sample, y=matNormalised, fill=cd133)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere normalised to no matrix") +
    ggtitle("Invasive ability of GIC clones buy CD133 status at day 3") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

invD7NormP = ggplot(data=invD7Norm[[1]], aes(x=sample, y=matNormalised, fill=cd133)) + 
    scale_fill_manual(values=c("yellow", "skyblue3")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere normalised to no matrix") +
    ggtitle("Invasive ability of GIC clones by CD133 status at day 7") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(invD3NormP, invD7NormP)

############################################## Plot the CD133 normlaised values #################################################

invD3NormCDP = ggplot(data=invD3Norm[[2]], aes(x=sample, y=cd133Norm, fill=sample)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", 'seagreen', 'plum2')) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to CD133 negative") +
    ggtitle("Invasive ability of GIC clones by CD133 status at day 3") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

invD7NormCDP = ggplot(data=invD7Norm[[2]], aes(x=sample, y=cd133Norm, fill=sample)) + 
    scale_fill_manual(values=c("slateblue", "maroon4", 'slategray1', 'palevioletred3')) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to CD133 negative") +
    ggtitle("Invasive ability of GIC clones by CD133 status at day 7") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(invD3NormCDP, invD7NormCDP)
####################################################################################################################################

invD3M = mean(invD3Norm[[2]]$cd133Norm)
invD3S = sd(invD3Norm[[2]]$cd133Norm)
invD7M = mean(invD7Norm[[2]]$cd133Norm)
invD7S = sd(invD7Norm[[2]]$cd133Norm)

invasionSummary = as.data.frame(rbind(c(invD3M, invD3S), c(invD7M, invD7S), c(1,0), c(1,0)))

invasionSummary$cd133 = c('pos', 'pos', 'neg', 'neg')
invasionSummary$day = c('day 3', 'day 7', 'day 3', 'day 7')
colnames(invasionSummary) = c('mean', 'sd', 'cd133', 'day')

invasionSummaryP = ggplot(data=invasionSummary, aes(x=day, y=mean, fill=cd133)) + 
    scale_fill_manual(values=c("aquamarine1", "darkgoldenrod1")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("") + ylab("Cell number relative to DMSO control") +
    ggtitle("Effect of temozolomide on subpopulations of GICs") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
invasionSummaryP
################################################# HAVE NOT MANUALLY CHECKED OUTLIERS YET ####################################################





############################################## Read in the ELDA results #################################################
# Clear memory
rm(list=ls())
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
setwd('~/Documents/Cell_biology/microscopy/ELDA/140417_elda_6clones/')
data = read.delim('140424_ELDA_output.txt')
data$Patient = c('030a', '030a', '039', '039', '035', '035', '020', '020', '041', '041', '034a')

eldaRawP = ggplot(data=data, aes(x=Group, y=Estimate, fill=Patient)) + 
    #scale_fill_manual(values=c("aquamarine1", "darkgoldenrod1")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(0.9)) +
    xlab("") + ylab("Sphere efficiency") +
    ggtitle("Raw ELDA data") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

percentData = as.data.frame(data[,2:4] ^ -1 *100)
percentData = cbind(percentData, data$Group)
percentData$Patient = c('030a', '030a', '039', '039', '035', '035', '020', '020', '041', '041', '034a')
colnames(percentData) = c('Lower', 'Estimate', 'Upper', 'Clone', 'Patient')

eldaPercent = ggplot(data=percentData, aes(x=Clone, y=Estimate, fill=Patient)) + 
    #scale_fill_manual(values=c("aquamarine1", "darkgoldenrod1")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(0.9)) +
    xlab("") + ylab("Percent sphere formation") +
    ggtitle("Sphere forming efficiency of GICs") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(eldaRawP, eldaPercent)