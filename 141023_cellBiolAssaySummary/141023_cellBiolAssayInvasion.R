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
# Divide by 1000 to make axis nicer
invasion$mean = rowMeans(invasion[,c(4:6)], na.rm=T) / 1000
invasion$sd = apply(invasion[,c(4:6)], 1, sd, na.rm=T) / 1000
colnames(invasion) = c("patient", "assayDate", "subpop", "rep1", "rep2", "rep3", "mean", "sd", "treatment")

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

spherePlotNorm = ggplot(dnNormSphere, aes(x=patient, y=normDN, fill=subpop)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=bw) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("PDGC") + ylab("Surface area relative to \nCD44-/CD133-") +
    ggtitle("Gliomasphere surface area at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
