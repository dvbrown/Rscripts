# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)

setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')

rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)

# Mung data into form for GSVA
rnaseqM = as.matrix(rnaseq)
sigs = list("CD133" = cd133Sig, "CD44" = cd44Sig)

result = gsva(rnaseqM, sigs, method="ssgsea", rnaseq=T, verbose=T)