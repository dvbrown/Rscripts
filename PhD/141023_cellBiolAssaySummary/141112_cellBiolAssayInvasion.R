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

setwd("~/Documents/Cell_biology/141023_summary/invasion/")
list.files()

# Intialise and write into database
db <- dbConnect(SQLite(), dbname="assaySummary.sqlite")
rawInvasion = read.delim("../141028_invasionSummary.txt")
invasion = read.delim("141112_mungInvasion.txt")

invasion$mean = rowMeans(invasion[,c(4:6)], na.rm=T)

bw = c("grey21", "grey82", "grey52", "grey97")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")

invPlot = ggplot(invasion, aes(x=patient, y=dnNorm, fill=subpop)) + 
    scale_fill_manual(values=color) + #guides(fill=FALSE) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("GSPC") + ylab("Invasive index relative to CD44-/CD133-") +
    ggtitle("Invasive index at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
invPlot

################## Invasion assay ############################
invSummary <- ddply(invasion, 'subpop', summarise,
                       N    = length(dnNorm), mean = mean(dnNorm),
                       geoMean = exp(mean(log(dnNorm))), geoSD = exp(sd(log(dnNorm))),
                       sd   = sd(dnNorm), se   = sd / sqrt(N) )

invSumPlot = ggplot(invSummary, aes(x=subpop, y=mean, fill=subpop)) + 
    #scale_fill_manual(values=color) + guides(fill=FALSE) +
    scale_fill_manual(values=bw) +  guides(fill=FALSE) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Invasive index relative to CD44-/CD133-") +
    ggtitle("Invasive index at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
invSumPlot

# pdf(file="./141112_invasionBw.pdf", useDingbats=F, height=12, width=18)
# invSumPlot
# dev.off()

anova(lm(mean ~ subpop + patient + treatment, data = rawInvasion)) # 0.0354165 *
TukeyHSD.aov(aov(mean ~ subpop + patient, data = rawInvasion), which="subpop") # CD44+/CD133-  CD44-/CD133- 0.0467581

anova(lm(dnNorm ~ subpop + patient, data = invasion)) # 0.04551 *
TukeyHSD.aov(aov(lm(dnNorm ~ subpop + patient, data = invasion)), which="subpop")

dbWriteTable(conn = db, name = "normalisedSphereArea", value = dnNormSphere, row.names = TRUE)
dbWriteTable(conn = db, name = "invasionSummary", value = invasionSummary, row.names = TRUE)
dbDisconnect(db)