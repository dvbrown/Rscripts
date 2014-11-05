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
setwd("~/Documents/Cell_biology/proliferation/Resazurin/141028_correct035/")
rawData = read.delim("pdgc035Reped.txt")
rawData = rawData[!is.na(rawData$subpop),]

backgroundMeanSD(rawData, 7)
