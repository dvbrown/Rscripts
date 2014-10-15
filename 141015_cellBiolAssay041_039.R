library(plyr)
source('140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

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
setwd("~/Documents/Cell_biology/microscopy/invasion/141014_clone041_039/")
rawData = read.delim("141015_InvasionRep.txt")

invasion = backgroundMeanSD(rawData)

# Plot Data
spherePlot = ggplot(data=invasion[invasion$treatment %in% FALSE,], aes(x=patient, y=mean, fill=subpop)) + 
    #scale_fill_manual(values=c("gold", "chartreuse4", "skyblue2", "forestgreen")) +
    scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Surface area of sphere without matrix") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
spherePlot

invasionNorm = normaliseMatrixCD133(invasion)
invasionSubpop = invasionNorm[[2]]

normPlot = ggplot(invasionSubpop, aes(x=patient, y=norm, fill=subpop)) + 
    #scale_fill_manual(values=c("skyblue2", "forestgreen")) +
    scale_fill_manual(values=c("lightgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Surface area realtive to double negative") +  # Set title
    # scale_x_continuous(breaks=pretty_breaks(5)) +
    scale_y_continuous(breaks = round(seq(0, 9, by = 1),1)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
normPlot

multiplot(spherePlot, normPlot)

# write.table(invasion, "./140927_004AND041Analysis/141008_invasionReplicates.txt", sep='\t')
# write.table(invasionNorm[[2]], "./140927_004AND041Analysis/141008_invasionRNormalised.txt", sep='\t')
# write.table(invasionNorm[[1]], "./140927_004AND041Analysis/141008_invasionRmatNormalised.txt", sep='\t')