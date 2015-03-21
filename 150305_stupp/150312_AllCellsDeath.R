# Analyse the dead cells
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Treatment']),]
    baseline = dataSorted[dataSorted$Treatment %in% baseline,]
    treatment = dataSorted[dataSorted$Treatment %in% value,]
    subtract =  treatment[,c(4)] - baseline[,c(4)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

plotPDGC = function(dataFrame, pdgc="MU035") {
    # This function takes a dataFrane of some data and returns a ggplot object specific to the cell line supplied as a string
    extract = dataFrame[,"PDGC"] %in% pdgc
    datPlot = dataFrame[extract,]
    cols = c("orange", "darkred", "darkblue", "forestgreen")
    p = ggplot(data=datPlot, aes(x=Treatment, y=value, fill=variable)) +
        geom_bar(stat="identity", position=position_dodge(), colour="black") + 
        ggtitle(paste("Viability after 7 days of treatment", pdgc, sep=" ")) +  scale_fill_manual(values=cols) + 
        xlab("Treatment") + ylab("Frequency") +
        theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return (p) }

setwd("~/Documents/Rscripts/150305_stupp/dat/")
list.files()

deadFreq = read.csv("150312_deadCells.csv", row.names=1)
deadLong = melt(deadFreq, id.vars=c("PDGC", "Treatment"))
# Subset only the dead cell column
dead = deadLong[deadLong$variable %in% "DeadCells",]

cols = c("darkgreen", "orange", "red", "yellow")
p1 = ggplot(data=dead, aes(x=PDGC, y=value, fill=Treatment)) +
        geom_bar(stat="identity", position=position_dodge(), colour="black") + 
        ggtitle("Viability after 7 days of treatment") +  scale_fill_manual(values=cols) + 
        xlab("Treatment") + ylab("Frequency") +
        theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Subtract the DMSO baseline
subtracted = rbind(subtractBaseline(dead, "dmso", "tmz"), subtractBaseline(dead, "dmso", "rad"),
                   subtractBaseline(dead, "dmso", "stupp"))

colours = c("orange", "red", "yellow")
p2 = ggplot(data=subtracted, aes(x=PDGC, y=subtract, fill=Treatment)) +
        geom_bar(stat="identity", position=position_dodge(), colour="black") + 
        ggtitle("Viability after 7 days of treatment") +  scale_fill_manual(values=colours) + 
        xlab("Treatment") + ylab("Frequency less DMSO") +
        theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(p1, p2, cols=1)

subtracted$subtract = 100 - subtracted$subtract
bioRep = ddply(subtracted[!subtracted$PDGC %in% "MU004",], .(Treatment), summarise, meanViability = mean(subtract, na.rm=T), 
               sdDiff = sd(subtract, na.rm=T), reps=length(subtract))
bioRep$seDiff = bioRep$sdDiff / (bioRep$reps)

ggplot(data=bioRep, aes(x=Treatment, y=meanViability, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("Viability after 7 days of treatment") +  scale_fill_manual(values=colours) + 
    xlab("Treatment") + ylab("Viability difference relative to DMSO") +
    scale_y_continuous(limits = c(0, 100)) +
    geom_errorbar(aes(ymin=meanViability-seDiff, ymax=meanViability+seDiff), width=.2, position=position_dodge(0.9)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
