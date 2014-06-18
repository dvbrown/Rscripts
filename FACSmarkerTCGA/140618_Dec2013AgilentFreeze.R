#try old Agilent dataSet 2013-12-10 freeze to check
library(GSVA)
library(coin)
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

statusClinical = function(clinicalDataFrame, vitalStatus) {
    status = clinicalDataFrame[,vitalStatus]
    clinicalDataFrame$status = ifelse(status == 'alive', 0, 1)
    return (clinicalDataFrame)
}

callMarkerSubtype <- function (signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # INCLUDE a double Negative
    signatureScore$subtype[signatureScore$CD133 < 0 & signatureScore$CD44 < 0] = "doubleNegative"
    signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}

setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
gem = read.delim('140110_agilentNoNulls.txt')
clinical = read.delim('../140109_clinicalDataTCGA.txt')
clinical = clinical[!is.na(clinical$patient.followups.followup.vitalstatus),]
clinical = statusClinical(clinical, 'patient.followups.followup.vitalstatus')

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))

# Extract the clinical data for the Agilent patients
matched = intersect(row.names(clinical), colnames(gem))
# Subset clinical data for intersect
gemM = gem[,matched]
clin = clinical[matched, ]

gemT = (as.matrix(gem))
gemT [is.na(gemT)] <- 0

# Call the subtypes with GSVA
# bigResult = gsva(gemT, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
# bigResult = t(bigResult$es.obs)

bigResult = gsva(gemT, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1, method='ssgsea')
bigResult = t(bigResult)

# Call a CD44 and CD133 subtype
cd133_44Sig = as.data.frame(bigResult[,c('CD133', 'CD44')])
verhaakSubtypeCall = callMarkerSubtype(cd133_44Sig, 0, 0)

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$patient.gender = as.factor(boundData$patient.gender)


#### Build survival object ####
data.surv = Surv(boundData$patient.daystodeath, event=boundData$status)

sur.fit = survfit(data.surv~subtype, boundData)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by Agilent',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue', 'orange'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44', 'double neg'),# 'Intermediate'), 
       col=c("red",'blue', 'orange'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype, subset=!boundData$subtype %in% "doubleNegative")
test
#text(locator(1),labels='p=0.0151', cex=1) #add the p-value to the graph