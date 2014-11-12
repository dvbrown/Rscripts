library(sqldf)
library(ggplot2)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd("~/Documents/Cell_biology/141023_summary/elda/")
list.files()

summariseByFactor = function (dataFrame, factor1) {
    require(plyr)
    # Takes a dataframe with factor information and computes summary statistics based on levels of as least 2 factors
    # factor 1 and 2 are characters
    result <- ddply(dataFrame, factor1, summarise,
                    N    = length(norm),
                    mean = mean(norm),
                    sd   = sd(norm),
                    se   = sd / sqrt(N) )
    return (result)
}

# Intialise and write into database
db <- dbConnect(SQLite(), dbname="assaySummary.sqlite")
elda = read.delim("141111_eldaSummary.txt")
bw = c("grey21", "grey82", "grey52", "grey97")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")

efficieny = elda[]
efficieny[,c(4:6)] = (efficieny[,c(4:6)] ^ -1) * 100

# Dividing by 1.96 gives SE
ci = efficieny[,6] - efficieny[,5]
efficieny$se = ci / 1.96
efficieny$sample = paste(efficieny$clone, efficieny$date, sep="_")
efficieny = efficieny[c(1:16),]

# Summarise data
eldaSummary <- ddply(efficieny, 'subpopulation', summarise,
                    N    = length(estimate), mean = mean(estimate),
                    geoMean = exp(mean(log(estimate))), geoSD = exp(sd(log(estimate))),
                    sd   = sd(estimate), se   = sd / sqrt(N) )




multiplot(eldaNormPlot, eldaPlot, eldaPlotRetain, eldaSummPlot, cols=2)