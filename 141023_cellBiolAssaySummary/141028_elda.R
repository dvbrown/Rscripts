library(sqldf)
library(ggplot2)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd("~/Documents/Cell_biology/141023_summary/")
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
summariseByFactor(efficieny, 'subpopulation')

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
    ggtitle("Double sorted cells") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=24))
eldaPlot

efficieny = efficieny[c(1:8),]
eldaPlotRetain = ggplot(efficieny, aes(x=clone, y=estimate, fill=subpopulation)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.2, position=position_dodge(0.9)) +
    xlab("GPSC line") + ylab("Sphere efficiency") +
    ggtitle("Double sorted cells") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
eldaPlotRetain

fit <- aov(estimate ~ subpopulation + clone, efficieny) 
summary(fit) # 0.248

# normalise to double negative
efficieny$norm = efficieny[,5] / 4.115226
efficieny$norm[5:8] = efficieny[c(5:8),5] / 6.142506

eldaNormPlot = ggplot(efficieny, aes(x=clone, y=norm, fill=subpopulation)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    xlab("GPSC line") + ylab("Sphere forming efficiency\n normalised to -/-") +
    ggtitle("Double sorted cells") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
eldaNormPlot

summData = summariseByFactor(efficieny, 'subpopulation')

eldaSummPlot = ggplot(summData, aes(x=subpopulation, y=mean, fill=subpopulation)) + 
    scale_fill_manual(values=color) +
    #scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("subpopulation") + ylab("Sphere forming efficiency") +
    ggtitle("Summary plot n = 2") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size=24))
eldaSummPlot

fit2 <- aov(norm ~ subpopulation, efficieny) 
summary(fit2) # 0.0964
TukeyHSD(fit2)

fit3 <- aov(estimate ~ subpopulation, efficieny) 
summary(fit3) # 0.27

multiplot(eldaNormPlot, eldaPlot, eldaPlotRetain, eldaSummPlot, cols=2)