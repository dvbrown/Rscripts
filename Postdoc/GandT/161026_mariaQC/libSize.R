library(plyr)
library(ggplot2)
library(reshape)
source('~/Code/Rscripts/Templates/multiplot.R')

inputDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/Counts_Maria/"
exportDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/"
setwd(inputDir)

samplesheet = read.csv("sampleSheetAll.csv", row.names = 1)
samplesheet = samplesheet[,c(1:4,6:8)]
rownames(samplesheet) = samplesheet[,5]

b1 = read.csv("Batch1/b1_htseq_counts_all.csv", row.names = 1)
b2 = read.csv("Batch2/b2_htseq_counts_all.csv", row.names = 1)
b3 = read.csv("Batch3/b3_htseq_counts_all.csv", row.names = 1)
b4 = read.csv("Batch4/b4_htseq_counts_all.csv", row.names = 1)

count = cbind(b1, b2, b3, b4)
rm(b1, b2, b3, b4)
write.table(count, "161026_geneExpMatrix.txt", sep="\t")

libSize = as.data.frame(colSums(count))
libMatrix = merge(count, libSize)