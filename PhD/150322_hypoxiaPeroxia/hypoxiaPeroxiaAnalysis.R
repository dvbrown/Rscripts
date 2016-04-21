# Analyse the dead cells
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Treatment']),]
    baseline = dataSorted[dataSorted$Treatment %in% baseline,]
    treatment = dataSorted[dataSorted$Treatment %in% value,]
    subtract =  treatment[,c(3:8)] - baseline[,c(3:8)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

setwd("~/Documents/Rscripts/150322_hypoxiaPeroxia/dat/")
list.files()
freq = read.csv("hypoxiaPeroxia.csv", row.names=1)
colnames(freq) = c("PDGC","Treatment","Alive","CD44+","CD44+/CD133+" ,"CD133+",   
                   "CD44-/CD133-","notDebris")
subtracted = read.csv("hypoxiaPeroxiaSubtracted.csv", row.names=1)
colnames(subtracted) = c("PDGC","Treatment","Alive","CD44+","CD44+/CD133+" ,"CD133+",   
                         "CD44-/CD133-","notDebris")