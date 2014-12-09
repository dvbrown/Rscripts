# Take the gene expression of RNA seq batch1 and attempt to get copy number

setwd("~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/copyNumber/")
cpm = read.delim("../GLMedgeR/131021_normalisedCPM.txt", row.names=1)
head(cpm)

# Load the genomic coordinates in a bed like format