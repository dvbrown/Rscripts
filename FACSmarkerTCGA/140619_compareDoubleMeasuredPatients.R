# Comapre the same patients that were RNAseq'd and Agilented
library(GSVA)
library(gplots)
library(RColorBrewer)
library(survival)
library(coin)
source("~/Documents/Rscripts/120704-sortDataFrame.R")
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
setwd('~/Documents/public-datasets/cancerBrowser/dedupAgilent/')

myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)
rnaSeq = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)
rnaseqGem = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)

agilent = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/clinical_dataDots.txt", row.names=1)
agilentGem = read.delim('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/genomicMatrix', row.names=1)

matched = intersect(colnames(rnaseqGem), colnames(agilentGem))

rnaseqGemM = rnaseqGem[,matched]
agilentGemM = agilentGem[,matched]

par(mfrow=c(2,1))
boxplot(agilentGemM, main='Matched Agilent')
boxplot(rnaseqGemM, main='Matched RNAseq')

#### The signature computation ####

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)

agilentM = as.matrix(agilentGemM)
rnaSeqM = as.matrix(rnaseqGemM)

# Change NAs to 0 as GSVA doesn't like it
agilentM[is.na(agilentM)] <- 0

#Fix the row.names of the clinical names to atch dfata
#row.names(agilentM) = paste(row.names(agilentM), '.01', sep='')

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), colnames(agilentM))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

# Call the subtypes with GSVA
resultAgilent = gsva(agilentM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
resultAgilent = t(resultAgilent$es.obs)
resultRNAseq = gsva(rnaSeqM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
resultRNAseq = t(resultRNAseq$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)

subtypeAgilent = bindGeneExprCIMPClinical(clin, resultAgilent, signatures)
subtypeAgilent = sort.dataframe(subtypeAgilent, 'colours')
write.table(subtypeAgilent, "output.txt", sep='\t')
subtypeAgilent = read.delim("output.txt", row.names=1)
subTypeHeatAgilent = as.matrix(subtypeAgilent[,signatures])

subtypeRNA = bindGeneExprCIMPClinical(clin, resultRNAseq, signatures)
subtypeRNA = sort.dataframe(subtypeRNA, 'colours')
write.table(subtypeRNA, "output.txt", sep='\t')
subtypeRNA = read.delim("output.txt", row.names=1)
subTypeHeatRNA = as.matrix(subtypeRNA[,signatures])

# Make heat map with Veerhaak subtype
par(mfrow=c(2,1))

heatmap.2(t(subTypeHeatAgilent), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=subTypeHeatAgilent$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeatAgilent), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

heatmap.2(t(subTypeHeatRNA), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=subTypeHeatRNA$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeatRNA), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")