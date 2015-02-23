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

####################### IO ###########################
setwd("~/Documents/Cell_biology/microenvironment/150216_doseResponse/")
list.files()
dasBos3 = read.delim("DasBos_day3.txt")
dasBos3$Group = paste(dasBos3$Patient, dasBos3$Drug, dasBos3$Conc, sep="_")
dasBos7 = read.delim("dasBos_day5.txt")
dasBos7$Group = paste(dasBos7$Patient, dasBos7$Drug, dasBos7$Conc, sep="_")
ruxIL3 = read.delim("RuxIL_day3.txt")
ruxIL5 = read.delim("ruxIL6.txt")