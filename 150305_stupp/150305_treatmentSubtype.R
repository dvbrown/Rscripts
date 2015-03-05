# A script to analyse the FACS experiment that I have done
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

setwd("~/Documents/facsData/flowJo/fortessa/150302_stupp/data/")
list.files()

dat = read.delim("150303_subpopTreatment.txt")
colnames(dat) = c("Sample","PDGC","Treatment","CD44+CD133-","CD44+CD133+","CD44-CD133+","CD44-CD133-")