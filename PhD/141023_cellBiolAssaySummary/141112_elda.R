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
elda = read.delim("141120_eldaSummary.txt")
bw = c("grey21", "grey82", "grey52", "grey97")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")

efficieny = elda[]
efficieny[,c(4:6)] = (efficieny[,c(4:6)] ^ -1) * 100

# Dividing by 1.96 gives SE
ci = efficieny[,6] - efficieny[,5]
efficieny$se = ci / 1.96
efficieny$sample = paste(efficieny$clone, efficieny$date, sep="_")
efficieny = efficieny[c(1:12),]

eldaPlot = ggplot(efficieny, aes(x=clone, y=estimate, fill=subpopulation)) + 
    scale_fill_manual(values=color) + #guides(fill=FALSE) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    #geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("GSPC") + ylab("Sphere forming efficiency") +
    ggtitle("Sphere forming efficiency at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
eldaPlot

# Summarise data
eldaSummary <- ddply(efficieny, 'subpopulation', summarise,
                    N    = length(estimate), mean = mean(estimate),
                    geoMean = exp(mean(log(estimate))), geoSD = exp(sd(log(estimate))),
                    sd   = sd(estimate), se   = sd / sqrt(N) )

eldaSumPlot = ggplot(eldaSummary, aes(x=subpopulation, y=mean, fill=subpopulation)) + 
    scale_fill_manual(values=color) + guides(fill=FALSE) +
    #scale_fill_manual(values=bw) +  guides(fill=FALSE) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Sphere forming efficiency") +
    ggtitle("Sphere forming efficiency at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
eldaSumPlot

# pdf(file="./141120_elda.pdf", useDingbats=F, height=12, width=18)
# eldaSumPlot
# dev.off()

anova(lm(estimate ~ subpopulation + clone, data = efficieny)) #  0.0304706 *
TukeyHSD.aov(aov(lm(estimate ~ subpopulation + clone, data = efficieny)), which="subpopulation")
# CD44+/CD133+-CD44-/CD133- 0.0261761