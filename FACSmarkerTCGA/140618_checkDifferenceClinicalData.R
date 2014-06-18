# compare survival parameters Agilent and RNA-seq

rnaSeq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)
rnaseqGem = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)

agilent = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/clinical_dataDots.txt", row.names=1)
agilentGem = read.delim('~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt', row.names=1)

matchedR = intersect(colnames(rnaseqGem), row.names(rnaSeq))
matchedA = intersect(row.names(agilentGem), row.names(agilent))

clinR = rnaSeq[matchedR, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender")]
clinA = agilent[matchedA, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender")]

par(mfrow=c(2,1))
hist(clinR$CDE_survival_time, breaks='FD', main='RNAseq', xlim=c(0,3000), ylim=c(0,50))
hist(clinA$CDE_survival_time, breaks='FD', main="Agilent", xlim=c(0,3000), ylim=c(0,50))
par(mfrow=c(1,1))

##### Analyse the set of patients measured by both RNAseq and Agilent ####