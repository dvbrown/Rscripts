# Classify the Agilent orginal samples from Broad firehose. I did some preprocessing here I forget what
library(GSVA)
library(gplots)
library(RColorBrewer)
library(survival)
library(coin)
source("~/Documents/Rscripts/120704-sortDataFrame.R")
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

data = read.delim('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/140110_agilentNoNulls.txt')
tcgaSigs = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)
cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)

agilentM = as.matrix(data)
# Change NAs to 0 as GSVA doesn't like it
agilentM[is.na(agilentM)] <- 0

#Fix the row.names of the clinical names to atch dfata
row.names(agilentM) = paste(row.names(agilentM), '.01', sep='')

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig))

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(agilentM))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

# Call the subtypes with GSVA
bigResult = gsva(agilentM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
bigResult = t(bigResult$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
verhaakSubtype = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
# write.table(verhaakSubtype, "./survival/140603_verhaakSubtypeAgilent.txt", sep='\t')

verhaakSubtype = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtype[,signatures])
myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")



verhaakSubtypeCall = callMarkerSubtype(verhaakSubtype, 0, 0)

############################################## Survival analysis #############################################
boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

############################################# Analysing the data for survival ##################################

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData)

par(mfrow=c(2,1))
plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by original Agilent',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype)#, subset=!boundData$subtype %in% "intermediate")
test


############################################# Make a double negative partition ##################################
callMarkerSubtype <- function (signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # Not having and intermediate case is also better for the Kaplan Myer curve
    signatureScore$subtype[signatureScore[,"CD133"] < 0 & signatureScore[,"CD44"] < 0] = "doubleNegative"
    signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}

boundData3Group = callMarkerSubtype(boundData,0,0)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData3Group$CDE_survival_time, event=boundData3Group$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData3Group)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by Agilent',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue', 'orange'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44', 'Double Negative'),# 'Intermediate'), 
       col=c("red",'blue', 'orange'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData3Group$subtype, subset=!boundData3Group$subtype %in% 'doubleNegative')
test

############################################## Don't set NAs to 0 in Agilent #############################################
agilentM = t(as.matrix(data))

# Remove patients with incomplete gene measurements
agilentMComplete = agilentM[complete.cases(agilentM),]
agilentMComplete = t(agilentMComplete)

#Fix the row.names of the clinical names to atch dfata
colnames(agilentMComplete) = paste(colnames(agilentMComplete), '.01', sep='')

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), colnames(agilentMComplete))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

# Call the subtypes with GSVA
bigResult = gsva(agilentMComplete, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
bigResult = t(bigResult$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
verhaakSubtype = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# The damn datatypes are not correct. Dump and read in object from file
write.table(verhaakSubtype, "output.txt", sep='\t')
# write.table(verhaakSubtype, "./survival/140603_verhaakSubtypeAgilent.txt", sep='\t')

verhaakSubtype = read.delim("output.txt", row.names=1)

subTypeHeat = as.matrix(verhaakSubtype[,signatures])

# Make heat map with Veerhaak subtype
heatmap.2(t(subTypeHeat), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype", 
          Colv=verhaakSubtype$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(verhaakSubtype$colours), labRow=colnames(subTypeHeat), xlab="Aglent samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")



verhaakSubtypeCall = callMarkerSubtype(verhaakSubtype, 0, 0)

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData)

par(mfrow=(c(2,1)))
plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by original Agilent',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype, subset=!boundData$subtype %in% "doubleNegative")
test
legend(locator(1), 'p-value = 0.034')

callMarkerSubtype <- function (signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # Not having and intermediate case is also better for the Kaplan Myer curve
    signatureScore$subtype[signatureScore[,"CD133"] < 0 & signatureScore[,"CD44"] < 0] = "doubleNegative"
    signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}
#### Try3 subtypes ####

verhaakSubtypeCall = callMarkerSubtype(verhaakSubtype, 0, 0)

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by original Agilent',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue', 'orange'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44', 'double negative'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype, subset=!boundData$subtype %in% "doubleNegative")
test
legend(locator(1), 'p-value (of CD44 vs CD133) = 0.11')
par(mfrow=(c(1,1)))
