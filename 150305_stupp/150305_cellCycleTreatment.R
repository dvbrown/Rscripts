# A script to to profile the cell cycles changes with treatment
require(ggplot2)
library(reshape)
source("./multiplot.R")

plotPDGC = function(dataFrame, pdgc) 
    {
    datPlot = dataFrame[dataFrame$pdgc %in% pdgc,]
    p = ggplot(data=ddatPlot, aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) + ggtitle(paste("PDGC", pdgc, sep=" ")) +  # Set title
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return (p)
}
p3 = plotPDGC(freqLong, "MU035")
p3

setwd("./150305_stupp/dat/")
col = c("darkgreen", "blue", "yellow", "darkred")
dat = read.delim("150308_cellCycleTreatment.txt")
View(dat)

freq = dat[,c(1:6)]
colnames(freq) = c("Sample", "PDGC", "Treatment", "G1", "G2", "S")
count = dat[,c(1:3,7:9)]

# Convert from wide to long
freqLong = melt(freq, id.vars=c("Sample", "PDGC", "Treatment"))