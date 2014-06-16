source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
setwd('~/Documents/Cell_biology/microscopy/invasion/140606_4clones/')

invD7 = read.delim('140616_processedOutputV2.txt')

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

# Plot Data
invD3StatsP = ggplot(data=invD3Stats[!invD3Stats$clone %in% c('030a', '034a'),], aes(x=sample, y=mean, fill=matrix)) + 
    scale_fill_manual(values=c("royalblue", "darkorange")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Invasive ability of GIC clones buy CD133 status at day 3") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

invD7StatsP = ggplot(data=invD7Stats[!invD7Stats$clone %in% c('030a', '034a'),], aes(x=sample, y=mean, fill=matrix)) + 
    scale_fill_manual(values=c("skyblue3", "yellow")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Invasive ability of GIC clones buy CD133 status at day 7") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(invD3StatsP, invD7StatsP)

############################################## Normalise for no matrix and CD133 #################################################

# Remove the unmatched recurrents
invD3Stats = invD3Stats[!invD3Stats$clone %in% c('030a', '034a', '20'),]
# invD7Stats = invD7Stats[!invD7Stats$clone %in% c('030a', '034a'),] For heatMap
invD7Stats = invD7Stats[!invD7Stats$clone %in% c('030a', '034a', '020'),]

invD3Norm = normaliseMatrixCD133(invD3Stats)
invD7Norm = normaliseMatrixCD133(invD7Stats)

# noMatrix = invD7Stats[invD7Stats$matrix %in%  FALSE,] This part is to calculate 020 neg for the heatmap
# noMatrix = noMatrix[c(1:7),]
# matrix = invD7Stats[invD7Stats$matrix %in% TRUE,]
# matrix$matNormalised = matrix$mean / noMatrix$mean
# matrix$matNormalisedSD = matrix$sd / noMatrix$sd

# Plot Data
invD3NormP = ggplot(data=invD3Norm[[1]], aes(x=sample, y=matNormalised, fill=clone)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", "forestgreen")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to no matrix") +
    ggtitle("Invasive ability of GIC clones \nby CD133 status at day 3") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

invD7NormP = ggplot(data=invD7Norm[[1]], aes(x=sample, y=matNormalised, fill=clone)) + 
    scale_fill_manual(values=c("yellow", "skyblue3", "darkseagreen2")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to no matrix") +
    ggtitle("Invasive ability of GIC clones \nby CD133 status at day 7") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# multiplot(invD3NormP, invD7NormP)

############################################## Plot the CD133 normlaised values #################################################

invD3NormCDP = ggplot(data=invD3Norm[[2]], aes(x=sample, y=cd133Norm, fill=sample)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", 'forestgreen')) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to CD133 negative") +
    ggtitle("Invasive ability of GIC clones \nby CD133 status at day 3") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

invD7NormCDP = ggplot(data=invD7Norm[[2]], aes(x=sample, y=cd133Norm, fill=sample)) + 
    scale_fill_manual(values=c("slateblue", "maroon4", 'darkseagreen2')) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to CD133 negative") +
    ggtitle("Invasive ability of GIC clones \nby CD133 status at day 7") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# multiplot(invD3NormCDP, invD7NormCDP)
multiplot(invD3NormP, invD7NormP, invD3NormCDP, invD7NormCDP, cols=2)
####################################################################################################################################

invD3M = mean(invD3Norm[[2]]$cd133Norm)
invD3S = sd(invD3Norm[[2]]$cd133Norm) / sqrt(3)
invD7M = mean(invD7Norm[[2]]$cd133Norm)
invD7S = sd(invD7Norm[[2]]$cd133Norm)/ sqrt(3)

invasionSummary = as.data.frame(rbind(c(invD3M, invD3S), c(invD7M, invD7S), c(1,0), c(1,0)))

invasionSummary$cd133 = c('positive', 'positive', 'negative', 'negative')
invasionSummary$day = c('day 3', 'day 7', 'day 3', 'day 7')
colnames(invasionSummary) = c('mean', 'sd', 'cd133', 'day')

invasionSummaryP = ggplot(data=invasionSummary, aes(x=day, y=mean, fill=cd133)) + 
    scale_fill_manual(values=c("forestgreen", "firebrick")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3, position=position_dodge(0.9)) +
    xlab("Date of measurement") + ylab("Surface area relative \nto no matrix control") +
    ggtitle("Invasive potential of CD133 \nsubpopulations of GICs") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
invasionSummaryP
