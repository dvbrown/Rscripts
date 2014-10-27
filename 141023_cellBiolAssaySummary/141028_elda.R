library(sqldf)
library(ggplot2)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd("~/Documents/Cell_biology/141023_summary/")
list.files()

# Intialise and write into database
db <- dbConnect(SQLite(), dbname="assaySummary.sqlite")
elda = read.delim("141023_eldaSummary.txt")
bw = c("grey21", "grey82", "grey52", "grey97")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")

efficieny = elda[]
efficieny[,c(4:6)] = (efficieny[,c(4:6)] ^ -1) * 100

# Dividing by 1.96 gives SE
ci = efficieny[,6] - efficieny[,5]
efficieny$se = ci / 1.96

efficieny$sample = paste(efficieny$clone, efficieny$date, sep="_")

eldaPlot = ggplot(efficieny, aes(x=sample, y=estimate, fill=subpopulation)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.2, position=position_dodge(0.9)) +
    xlab("GPSC line") + ylab("Sphere forming efficiency") +
    ggtitle("GPSC 041") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=24))
eldaPlot