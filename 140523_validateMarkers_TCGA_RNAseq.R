# Classify the RNA-seq TCGA samples provide they are different patients than Agilent

setwd('~/Documents/public-datasets/cancerBrowser/')

agilent = read.delim('./TCGA_GBM_G4502A_07_2-2014-05-02/genomicMatrix')
rnaSeq = read.delim('./TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix')