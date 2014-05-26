# Classify the RNA-seq TCGA samples provide they are different patients than Agilent

setwd('~/Documents/public-datasets/cancerBrowser/')

# Remove from Agilent data the RNAseq cases
agilent = read.delim('./TCGA_GBM_G4502A_07_2-2014-05-02/genomicMatrix', row.names=1)
rnaSeq = read.delim('./TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix', row.names=1)

# Measure the patient overlap between Agilent and HiSeq
agilentNames = colnames(agilent)
# rm(agilent)
rnaSeqNames = colnames(rnaSeq)
head(agilentNames)
head(rnaSeqNames)
duplicatedPatients = intersect(agilentNames, rnaSeqNames)
# The number of unique genes
168 - length(intersect(agilentNames, rnaSeqNames))

