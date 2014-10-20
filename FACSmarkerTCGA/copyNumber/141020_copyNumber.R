library(sqldf)

setwd("~/Documents/public-datasets/TCGA/copyNumber//")
list.files()

anno = read.delim("TCGA_GBM_PCT_CNA_Annotation.txt", row.names=1)
triplet = read.delim("TCGA_GBM_PCT_CNA_Physical_triplet.txt", row.names=1)
summary = read.delim("TCGA_GBM_PCT_CNA_Summary.txt", row.names=1)