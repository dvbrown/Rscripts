# Classify the RNA-seq TCGA samples provide they are different patients than Agilent

setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/')

rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)
