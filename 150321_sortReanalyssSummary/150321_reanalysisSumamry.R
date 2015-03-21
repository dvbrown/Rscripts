# Thge script wll summarise the 2 experiments I did that sorted and then reanalysed PDGCs
library(ggplot2)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Subpopulation']),]
    baseline = dataSorted[dataSorted$Subpopulation %in% baseline,]
    treatment = dataSorted[dataSorted$Subpopulation %in% value,]
    subtract =  treatment[,c(4:7)] - baseline[,c(4:7)]
    result = cbind(treatment[,c(1:3)], subtract)
    return (result)
}

plotPDGC = function(dataFrame, pdgc="MU035") {
    # This function takes a dataFrane of some data and returns a ggplot object specific to the cell line supplied as a string
    extract = dataFrame[,"PDGC"] %in% pdgc
    datPlot = dataFrame[extract,]
    cols = c("orange", "darkred", "darkblue", "forestgreen")
    p = ggplot(data=datPlot, aes(x=Subpopulation, y=value, fill=variable)) +
        geom_bar(stat="identity", position=position_dodge(), colour="black") + 
        ggtitle(paste("Sorting an reanalysis after 7 days", pdgc, sep=" ")) +  scale_fill_manual(values=cols) + 
        xlab("Subpopulation") + ylab("Percent difference \nrelative to mixed population") +
        theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return (p) }

setwd("~/Documents/Rscripts/150321_sortReanalyssSummary/")
freq = read.csv("rawData.csv", row.names=1)

freqBase = rbind(subtractBaseline(freq, "mixed", "doubleNeg"), subtractBaseline(freq, "mixed", "CD44"),
                 subtractBaseline(freq, "mixed", "CD133"), subtractBaseline(freq, "mixed", "doublePos"))

freqBase = freqBase[order(freqBase[,'PDGC'], freqBase[,'Subpopulation']),]