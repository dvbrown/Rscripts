# Note that in this experiment I beleive that clone #035 has been swapped with #039
# See my lab book #2 pg 95

source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

############################################## Analyse invasion assay ###############################################
backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data take the mean and sd
    dataFrame$mean = rowMeans(dataFrame[,c(4:6)], na.rm=T)
    dataFrame$sd = apply(dataFrame[,c(4:6)], 1, sd, na.rm=T)
    return (dataFrame)
}

normaliseMatrixDN = function(dataFrame) {
    # Normalise Matrix first
    noMatrix = dataFrame[dataFrame$treatment %in%  FALSE,]
    matrix = dataFrame[dataFrame$treatment %in% TRUE,]
    matrix$matNormalised = matrix$mean / noMatrix$mean
    matrix$matNormalisedSD = matrix$sd / noMatrix$sd
    # Normalise doubleNeg
#     negative = matrix[matrix$subpop %in% 'CD44-/CD133-',]
#     theRest = matrix[!matrix$subpop %in% 'CD44-/CD133-',]
#     theRest$norm = theRest$mean / negative$mean
#     # Return both dataframes
#     result = list(matrix, theRest)
    return (matrix)
}

setwd("~/Documents/Cell_biology/microscopy/invasion/141014_clone041_039/")
rawData = read.delim("141015_InvasionRep.txt")
invasion = backgroundMeanSD(rawData)

bw = c("grey21", "grey82", "grey52", "grey97")
colour = c("gold", "chartreuse4", "skyblue2", "forestgreen")

# Plot Data
spherePlot = ggplot(data=invasion[invasion$treatment %in% FALSE,], aes(x=patient, y=mean, fill=subpop)) + 
    #scale_fill_manual(values=colour) +
    scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Surface area of sphere without matrix") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
spherePlot

invasionDN = normaliseMatrixDN(invasion)

# Have to manually make the normalisation
invasionDN$norm = 1
invasionDN[2,14] = invasionDN[2,12] / invasionDN[1,12]
invasionDN[4,14] = invasionDN[4,12] / invasionDN[3,12]
invasionDN[5,14] = invasionDN[5,12] / invasionDN[3,12]
invasionDN

normPlot = ggplot(invasionDN, aes(x=patient, y=norm, fill=subpop)) + 
    #scale_fill_manual(values=colour) +
    scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Surface area realtive to double negative") +  # Set title
    # scale_x_continuous(breaks=pretty_breaks(5)) +
    scale_y_continuous(breaks = round(seq(0, 9, by = 1),1)) +
    #geom_hline(yintercept=1) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
normPlot

multiplot(spherePlot, normPlot)

write.table(invasion, "./141015_meanSDInvasion.txt", sep='\t')
write.table(invasionDN, "./141015_invasionRNormalised.txt", sep='\t')

#rm(list=ls())

############################################## Analyse growth assay ###############################################
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

setwd("~/Documents/Cell_biology/proliferation/Resazurin/141007_3clones/")
rawData = read.delim("141015_growthDayRepBlanked.txt")

summarisedData = backgroundMeanSD(rawData)
bw = c("grey21", "grey82", "grey52", "grey97")
colour = c("chartreuse4", "gold", "skyblue2", "orangered1")

# Plot the raw results
growthPlot7 = ggplot(summarisedData[summarisedData$treatment %in% 'DMSO',], 
                     aes(x=patient, y=mean, fill=subpop)) + 
    #scale_fill_manual(values=colour) +
    scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("PDAC") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

# Analyse the temozolomide
temo = calcDMSOcontrol(summarisedData)

tmzPlot7 = ggplot(data=temo, aes(x=patient, y=dmsoCorrected, fill=subpop)) + 
    #scale_fill_manual(values=colour) +
    scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("PDAC") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Temozolomide sensitivty at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

# multiplot(growthPlot7, tmzPlot7)

# write.table(temo, './141015_processedData.txt', sep='\t')

############################################## Analyse ELDA assay ###############################################
#rm(list=ls())
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
setwd("~/Documents/Cell_biology/microscopy/ELDA/141007/")
rawData = read.delim("141015_eldaOutData.txt")
rawData = rawData[c(1:8),]

bw = c("grey21", "grey82", "grey52", "grey97")
colour = c("chartreuse4", "gold", "skyblue2", "orangered1")

# Transform to efficiencey
transformedData = rawData[]
transformedData[,c(2:4)] = (rawData[,c(2:4)] ^ -1) * 100

# Dividing by 1.96 gives SE
ci = transformedData[,3] - transformedData[,2]
transformedData$se = ci / 1.96
transformedData$sd = sqrt(96) * transformedData$se

eldaPlot = ggplot(transformedData, aes(x=patient, y=Estimate, fill=subpopulation)) + 
    #scale_fill_manual(values=colour) +
    scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=Estimate-se, ymax=Estimate+se), width=.2, position=position_dodge(0.9)) +
    xlab("GPSC line") + ylab("Sphere forming efficiency") +
    ggtitle("GPSC 041") + scale_y_continuous(breaks = round(seq(0, 100, by = 5),2)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

svg("~/Documents/Cell_biology/proliferation/Resazurin//141007_3clones/141015_allExptsTogether.svg", width=11.69, height=8.27)
multiplot(growthPlot7, tmzPlot7, eldaPlot, spherePlot, normPlot, cols=2)
dev.off()