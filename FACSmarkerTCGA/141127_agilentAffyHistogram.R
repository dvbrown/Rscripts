# A script to demonstrate Agilent is better than Affymetrix

#Load Agilent
setwd('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/')
agilent = read.delim("genomicMatrix", row.names=1)

# Load Affy
affy = read.delim('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_u133a-2014-08-22/genomicMatrix', row.names=1)

patients = intersect(colnames(affy), colnames(agilent))
# Take patients in common
agilent = agilent[,patients]
affy = affy[,patients]

affyGAPDH = as.data.frame(t(affy["GAPDH",]))
affyGAPDH$gene = "GAPDH"
affyBact = as.data.frame(t(affy["ACTB",]))
affyBact$gene = "ACTB"
colnames(affyGAPDH) = c("expression", "gene")
colnames(affyBact) = c("expression", "gene")
affyGene = rbind(affyGAPDH, affyBact)


agilentGAPDH = as.data.frame(t(agilent["GAPDH",]))
agilentGAPDH$origin = "Agilent"
