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

##### Box plot the data frpm RNAseq and microarray ####

par(mfrow=c(2,1))
boxplot(t(agilentGem[c(1:133)]), main="Agilent", col=rainbow(133))
boxplot(rnaseqGem[c(1:133)], main="RNAseq", col=rainbow(133))
par(mfrow=c(1,1))

# Read in the agilent file form Firehose
originalAgilent = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM//20131210_dataReformatting//dataRearranging/140109_agilent.txt')
originalAgilentM = originalAgilent[,matchedA]
data = read.delim('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/140110_agilentNoNulls.txt')

par(mfrow=c(2,1))
boxplot(originalAgilent[c(1:133)], main="Original Agilent from Firehose", col=rainbow(133))
boxplot(data[c(1:133)], main="Original Agilent from Firehose no Nulls", col=rainbow(133))
par(mfrow=c(1,1))


##### Analyse the set of patients measured by both RNAseq and Agilent ####