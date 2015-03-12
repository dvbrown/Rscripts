# Analyse the dead cells
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Treatment']),]
    baseline = dataSorted[dataSorted$Treatment %in% baseline,]
    treatment = dataSorted[dataSorted$Treatment %in% value,]
    subtract =  treatment[,c(3:7)] - baseline[,c(3:7)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

setwd("~/Documents/Rscripts/150305_stupp/dat/")
list.files()
col = c("darkgreen", "orange", "red", "yellow")
deadFreq = read.csv("150312_deadCells.csv", row.names=1)
deadFreq