# Classify the RNA-seq TCGA samples provide they are different patients than Agilent
library(GSVA)
library(gplots)

setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')

rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)

# Mung data into form for GSVA
rnaseqM = as.matrix(rnaseq)
sigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

result = gsva(rnaseqM, sigs, method="ssgsea", rnaseq=T, verbose=T)
result = t(result)

par(mfrow=c(2,1))
hist(result[,1], breaks='FD', main="CD133", xlim=c(0,1))
hist(result[,2], breaks='FD', main="CD44", xlim=c(0,1))
par(mfrow=c(1,1))

heatmap.2(result, Colv=NA, cexRow=0.5, cexCol=0.9, main="ssGSEA FACS markers", scale="column", keysize=1, trace="none")
heatmap.2(result, Colv=NA, cexRow=0.5, cexCol=0.9, main="ssGSEA FACS markers", keysize=1, trace="none")