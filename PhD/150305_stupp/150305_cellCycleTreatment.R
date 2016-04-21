# A script to to profile the cell cycles changes with treatment
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

extract = freqLong[,"PDGC"] %in% "MU035"

plotPDGC = function(dataFrame, pdgc) {
    extract = dataFrame[,"PDGC"] %in% pdgc
    datPlot = dataFrame[extract,]
    p = ggplot(data=datPlot, aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=cols) + ggtitle(paste("PDGC", pdgc, sep=" ")) +  # Set title
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return (p) }

setwd("./150305_stupp/dat/")
cols = c("darkgreen", "orange", "red", "yellow")
dat = read.delim("150308_cellCycleTreatment.txt")
View(dat)

freq = dat[,c(1:6)]
colnames(freq) = c("Sample", "PDGC", "Treatment", "G1", "G2", "S")
count = dat[,c(1:3,7:9)]

# Convert from wide to long
freqLong = melt(freq, id.vars=c("Sample", "PDGC", "Treatment"))

# Make the plots
p35 = plotPDGC(freqLong, "MU035")
p20 = plotPDGC(freqLong, "MU020")
p4 = plotPDGC(freqLong, "MU004")
p39 = plotPDGC(freqLong, "MU039")

multiplot(p4, p20, p35, p39, cols=2)