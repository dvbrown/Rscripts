library(sqldf)
library(ggplot2)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

summariseByFactor <- function (dataFrame, factor1, factor2, factor3) {
    require(plyr)
    # Takes a dataframe with factor information and computes summary statistics based on levels of as least 2 factors
    # factor 1 and 2 are characters
    result <- ddply(dataFrame, c(factor1, factor2, factor3), summarise,
                    N    = length(mean),
                    mean = mean(mean),
                    sd   = sd(mean),
                    se   = sd / sqrt(N) )
    return (result)
}

combineExperiments <- function (dataFrame) {
  # Run consective normalisations then bind and trim the result
  four = summariseByFactor(dataFrame, "patient", "assayDate", "treatment")
#   twenty = summaryFunction(dataFrame, "#020", assayDate)
#   twenty8 = summaryFunction(dataFrame, "#028")
#   thirty5 = summaryFunction(dataFrame, "#035")
#   thrity9 = summaryFunction(dataFrame, "#039")
#   # 041 was done twice, need to separate
#   fourty1 = summaryFunction(dataFrame[c(1:34),], "#041")
#   fourty2 = summaryFunction(dataFrame[c(35:38),], "#041")
#   rawData = rbind(four, twenty, thirty5, thrity9, fourty1, fourty2)
#   # Get rid of the duplicate #041 readings
#   deDup = rawData[c(1:15,17),]
  return (four)
}
combineExperiments(invasion)

setwd("~/Documents/Cell_biology/141023_summary/")
list.files()

# Intialise and write into database
db <- dbConnect(SQLite(), dbname="assaySummary.sqlite")
invasion = read.delim("141023_invasionSumRound.txt")
invasion$mean = rowMeans(invasion[,c(4:6)], na.rm=T)
invasion$sd = apply(invasion[,c(4:6)], 1, sd, na.rm=T)
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
    xlab("PDGC") + ylab("Fluorescent intensity") +
    ggtitle("Growth at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
