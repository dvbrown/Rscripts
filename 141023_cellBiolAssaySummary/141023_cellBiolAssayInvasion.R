library(sqldf)
library(ggplot2)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

calcInvNormalised = function(dataFrame, patientName) {
    # Normalises a patient by the double negative subpopulation set to 1
    # Patient is a charaacter string of patient eg #035
    # Extract all cases of the individual patient
    patient = dataFrame[dataFrame[,"patient"] %in% patientName,]
    patient = patient[patient[,"treatment"] %in% "FALSE",]
    # Extract the double negative
    dn = patient[patient[,"subpop"] %in% "CD44-/CD133-",]
    otherSample = patient
    otherSample$normDN = otherSample$mean / dn$mean
    return (otherSample)
}

combineExperiments <- function (dataFrame, summaryFunction) {
  # Run consective normalisations then bind and trim the result
  four = summaryFunction(dataFrame, "#004")
  twenty = summaryFunction(dataFrame, "#020")
  #thirty5 = summaryFunction(dataFrame, "#035")
  thrity9 = summaryFunction(dataFrame, "#039")
  rawData = rbind(four, twenty, thrity9)
  # Get rid of the duplicate #041 readings
  return (rawData)
}

setwd("~/Documents/Cell_biology/141023_summary/")
list.files()

# Intialise and write into database
db <- dbConnect(SQLite(), dbname="assaySummary.sqlite")
invasion = read.delim("141028_invasionSummary.txt")
invasion = invasion[c(1:26,28:31),]
# Divide by 1000 to make axis nicer
invasion$mean = rowMeans(invasion[,c(4:6)], na.rm=T) / 1000
invasion$sd = apply(invasion[,c(4:6)], 1, sd, na.rm=T) / 1000

bw = c("grey21", "grey82", "grey52", "grey97")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")

################## Invasion assay ############################
# Plot the raw results
spherePlot = ggplot(invasion[invasion$treatment %in% 'FALSE',], 
                    aes(x=patient, y=mean, fill=subpop)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("PDGC") + ylab("Surface area (uM)") +
    ggtitle("Gliomasphere surface area at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

# Normalise by double negative. As #035 has no DN must discard it
dnNormSphere = combineExperiments(invasion, calcInvNormalised)
# Need to mung #041 myself
fourty = read.delim("141023_invasion041.txt")
dnNormSphere = rbind(dnNormSphere, fourty[fourty$treatment %in% 'FALSE',])

spherePlotNorm = ggplot(dnNormSphere, aes(x=patient, y=normDN, fill=subpop)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("PDGC") + ylab("Surface area relative to \nCD44-/CD133-") +
    ggtitle("Gliomasphere surface area at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

################## matrix ############################
normaliseNoMat= function(dataFrame) {
    matrixFalse = dataFrame[dataFrame[,"treatment"] %in% "FALSE",]
    matrixTrue = dataFrame[dataFrame[,"treatment"] %in% "TRUE",]
    matrixTrue$noMatNorm = matrixTrue$mean / matrixFalse$mean
    return (matrixTrue)
}

summariseByFactor = function (dataFrame, factor1, factor2) {
    require(plyr)
    # Takes a dataframe with factor information and computes summary statistics based on levels of as least 2 factors
    # factor 1 and 2 are characters
    result <- ddply(dataFrame, c(factor1, factor2), summarise,
                    N = length(dnNorm), mean = mean(dnNorm), sd   = sd(dnNorm), se   = sd / sqrt(N) )
    return (result)
}

matNorm = normaliseNoMat(invasion)
# write.table(matNorm, './invasion/141028_mungInvasion.txt', sep='\t')
matNorm = read.delim("./invasion/141028_mungInvasion.txt")

matrixNorm= ggplot(matNorm, aes(x=patient, y=dnNorm, fill=subpop)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("PDGC") + ylab("Invasion relative to \nCD44-/CD133-") +
    ggtitle("Gliomasphere invasion at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

invasionSummary = summariseByFactor(matNorm, 'subpop', 'treatment')

invSumPlot = ggplot(invasionSummary, aes(x=subpop, y=mean, fill=subpop)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("PDGC") + ylab("Invasion relative to \nCD44-/CD133-") +
    ggtitle("Gliomasphere invasion at day 7") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=24))
invSumPlot

multiplot(spherePlot, spherePlotNorm, matrixNorm ,invSumPlot, cols=2)