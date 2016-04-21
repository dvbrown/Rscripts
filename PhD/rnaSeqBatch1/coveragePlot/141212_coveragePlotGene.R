# A script that draws a coverage plot of a selected gene

library(ggbio)
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(cummeRbund)
setwd('~/Documents/RNAdata/RNAseqAnalysis/121116_cuffOutFull/')
library(BSgenome.Hsapiens.UCSC.hg19)
#build the gene model plots for your gene of choice
txdb=makeTranscriptDbFromUCSC(genome='hg19',tablename='ensGene')