source('140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

############################################## IO ###############################################
setwd("~/Documents/Cell_biology/proliferation/Resazurin/")
pat004Data = read.delim("140821_004/pdac004Day7Label.txt", row.names=1)