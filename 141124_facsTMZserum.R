library(ggplot2)
library(plyr)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd("~/Documents/Cell_biology/proliferation/Resazurin/141124_tmzFcsFACS/")
list.files()

growth = read.delim("141124_tmzFcsFACS.txt")
growth$sample = paste(growth[,"patient"], growth[,"treatment"], sep="_")

# Subtract background, take mean and sd
growth[,c(4:6)] = growth[,c(4:6)] - rowMeans(growth[17,c(4:6)])
growth$mean = rowMeans(growth[,c(4:6)], na.rm=T)
growth$sd = apply(growth[,c(4:6)], 1, sd, na.rm=T)
growth$cv = growth$sd / growth$mean * 100
write.table(growth, "141125_facsTMZserum.txt", sep='\t')