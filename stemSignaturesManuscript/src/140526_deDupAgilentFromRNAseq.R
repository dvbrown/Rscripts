agilent = read.delim('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/genomicMatrix', row.names=1)
rnaSeq = read.delim('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix', row.names=1)

# Measure the patient overlap between Agilent and HiSeq
agilentNames = colnames(agilent)
rnaSeqNames = colnames(rnaSeq)
duplicatedPatients = intersect(agilentNames, rnaSeqNames)
uniquePatients = setdiff(agilentNames, rnaSeqNames)
# Subset out the duplicates
agilentDeDup = t(agilent)
# Return the unique cases to Agilent
agilentDeDup1 = agilentDeDup[uniquePatients,]
write.table(agilentDeDup1, file="~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", sep="\t")