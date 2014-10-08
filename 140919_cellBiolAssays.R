library(plyr)
source('140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

summaryStats = function(dataFrame, groups) {
    # The dataFrame with the readings
    # a vector of strings with the column names you wish to group by
    # Measurement is the numerical numbers you're interested in
    result = ddply(dataFrame, groups, summarise,
                   N    = sum(!is.na(fluorescence)),
                   mean = mean(fluorescence, na.rm=TRUE),
                   sd   = sd(fluorescence),
                   se   = sd / sqrt(N) )
    return (result)
}

############################################## IO ###############################################
setwd("~/Documents/Cell_biology/proliferation/Resazurin/")
pat004Data = read.delim("140821_004/pdac004Day7Label.txt", row.names=1)
pat041Data = read.delim("140919_gpsc041/pdac041Day7Label.txt", row.names=1)

rawData = rbind(pat004Data, pat041Data)
# Subtract background and take mean and SD
growthData = summaryStats(rawData, c('patient', 'subpop', 'treatment'))

#### Plot the raw results ####
growthPlot7 = ggplot(growthData[growthData$treatment %in% 'DMSO',], 
                     aes(x=patient, y=mean, fill=subpop)) + 
    #scale_fill_manual(values=c("darkorange", "royalblue", "forestgreen", "red")) +
    scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("PDAC") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))


############################################## Plot DMSO controlled  scores ###############################################
# Remove 041 DMSO as the TMZ was contaminated
growthData = growthData[c(1:10),]
temo = calcDMSOcontrol(growthData)

tmzPlot7 = ggplot(data=temo, aes(x=patient, y=dmsoCorrected, fill=subpop)) + 
    #scale_fill_manual(values=c("gold", "chartreuse4", "skyblue2", "forestgreen")) +
    scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("PDAC") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Temozolomide sensitivty at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

multiplot(growthPlot7, tmzPlot7)

# write.table(growthData, './140919_gpsc041/141008_summarisedGrowth.txt', sep='\t')
# write.table(temo, './140919_gpsc041/141008_summarisedTMZ.txt', sep='\t')

############################################## Analyse invasion assay ###############################################
rm(list=ls())
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data take the mean and sd
    dataFrame$mean = rowMeans(dataFrame[,c(4:6)], na.rm=T)
    dataFrame$sd = apply(dataFrame[,c(4:6)], 1, sd, na.rm=T)
    return (dataFrame)
}

normaliseMatrixCD133 = function(dataFrame) {
    # Normalise Matrix first
    noMatrix = dataFrame[dataFrame$treatment %in%  FALSE,]
    matrix = dataFrame[dataFrame$treatment %in% TRUE,]
    matrix$matNormalised = matrix$mean / noMatrix$mean
    matrix$matNormalisedSD = matrix$sd / noMatrix$sd
    # Normalise doubleNeg
    negative = matrix[matrix$subpop %in% 'CD44-/CD133-',]
    theRest = matrix[!matrix$subpop %in% 'CD44-/CD133-',]
    theRest$norm = theRest$mean / negative$mean
    # Return both dataframes
    result = list(matrix, theRest)
    return (result)
}

############################################## IO ###############################################
setwd("~/Documents/Cell_biology/microscopy/invasion/")
pat004Data = read.delim("140905_clone004/141001_004Rep.txt")
pat041Data = read.delim("140926_invasion041/141001_output041Rep.txt")

rawData = rbind(pat004Data, pat041Data)
# Remove the  #041 CD44+/CD133- case that has no control
invasion = backgroundMeanSD(rawData)
invasion = invasion[c(1:10),]

# Plot Data
spherePlot = ggplot(data=invasion[invasion$treatment %in% FALSE,], aes(x=patient, y=mean, fill=subpop)) + 
    scale_fill_manual(values=c("gold", "chartreuse4", "skyblue2", "forestgreen")) +
    #scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Surface area of sphere without matrix") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
spherePlot

invasionNorm = normaliseMatrixCD133(invasion)
invasionSubpop = invasionNorm[[2]]

normPlot = ggplot(invasionSubpop, aes(x=patient, y=norm, fill=subpop)) + 
    scale_fill_manual(values=c("skyblue2", "forestgreen")) +
    #scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Surface area of sphere without matrix") +  # Set title
    # scale_x_continuous(breaks=pretty_breaks(5)) +
    scale_x_continuous(breaks = round(seq(min(invasionSubpop$norm), max(invasionSubpop$norm), by = 0.5),1)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
normPlot