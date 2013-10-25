setwd('~/Documents/RNAdata/')

data = read.delim('RNAseqProgress.txt')
da = data[c(1:7),c(4,5,6,7,8,9,10,11,12,14)]
da[4,3] = as.factor('primary')
da = da[c(1,3,4,5,6,7),]
colnames(da) = c('Clone', 'Passage', 'Origin', 'Group', 'Survival', 'Age', 'cellNo','CD133', 'RNA', 'cDNA')

da$rna = as.integer(da$RNA)
da$cellNo = as.integer(da$cellNo)

colors = c('darkgreen', 'lightgreen','springgreen','cyan', 'blue1', 'lightblue')

par(mfrow=c(2,2), las=2, cex=1.25)

barplot(da$Survival, main='Patient Survival', ylab='Survival (months)', names.arg=da$Clone, xlab='Patient ID',
        col=colors)

barplot(da$Age, main='Patient Age', ylab='Age (years)', names.arg=da$Clone, xlab='Patient ID',
        col=colors)

barplot(da$Passage, main='Passage number', ylab='Passage number', names.arg=da$Clone, xlab='Patient ID',
        col=colors, ylim=c(0,10))

barplot(da$CD133, main='Patient CD133+', ylab='CD133 positive %', names.arg=da$Clone, xlab='Patient ID',
        col=colors, ylim=c(0,100))

######### New plots ##############
par(mfrow=c(2,2), cex=1.25)
barplot(as.integer(da$cellNo), main='Cell number', ylab="Cell number (000's)", names.arg=da$Clone, xlab='Patient ID',
        col=colors)

barplot(da$RNA, main='RNA extracted', ylab='Concentration (ng)', names.arg=da$Clone, xlab='Patient ID',
        col=colors)

barplot(da$cDNA, main='Library yield', ylab='Concentration (ng)', names.arg=da$Clone, xlab='Patient ID',
        col=colors, ylim=c(0,45))

################################## RNAseQC plots ##################################
rnaSEqc = read.delim('danBatch1/rnaSeQC_reports/131008_rnaseQCSummary.txt', row.names=1)
rnaSEqc = rnaSEqc[,c(1,3,4,6,7,10,11,12,13,14)]
par(mfrow=c(2,2), cex=1.1)

barplot(rnaSEqc$TotalReadsSequenced, main='Total reads sequenced', ylab='Reads (millions)',
        col=colors, names.arg=row.names(rnaSEqc), las=2, cex.names=0.8)

barplot(rnaSEqc$Mapped_Unique, main='Uniquely mapped reads', ylab='Reads (millions)',
        col=colors, names.arg=row.names(rnaSEqc), las=2, ylim=c(0,70), cex.names=0.8)

barplot(rnaSEqc$Duplication_Rate_, main='PCR and optical duplicates', ylab='Duplication rate',
        col=colors, names.arg=row.names(rnaSEqc), las=2, ylim=c(0,1), cex.names=0.8)

barplot(rnaSEqc$Intragenic_Rate, main='Reads located in coding reigons', ylab='Intragenic rate',
        col=colors, names.arg=row.names(rnaSEqc), las=2, ylim=c(0,1), cex.names=0.8)
