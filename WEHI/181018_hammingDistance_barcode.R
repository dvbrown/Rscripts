setwd('/Users/brown.d/Desktop/')
library(DNABarcodes)

dat = read.delim("barcodes.txt",header = F)
vec = as.vector(dat$V2)
names(vec) = dat$V1
vec

analyse.barcodes(vec)