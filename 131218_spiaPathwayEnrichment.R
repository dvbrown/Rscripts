# Measuring pathway enrichment using SPIA in RNA-seq batch1
library(SPIA)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/')

# Use the liberal differentially expressed genes
deGenes = read.delim('131021_shortVSlongLiberalDE.txt')