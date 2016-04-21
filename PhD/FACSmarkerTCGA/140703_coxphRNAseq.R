library(survival)
library(coin)
library(ggplot2)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

############################################## IO and general munging #############################################

verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140530_liberalSignatureScores2SD.txt", row.names=1)

verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS")]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

# Call the FACS subtype
verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", 'CDE_chemo_adjuvant_tmz', 'CDE_chemo_tmz')]

############################################## bind the clinical and subtyping info together #############################################

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)

############################################# Analysing the data for survival ##################################
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)
xtabs(~ CDE_DxAge + subtype, data=boundData)

coxObj = coxph(data.surv ~ CDE_DxAge + subtype, data=boundData, na.action="na.omit")
summary(coxObj)

# Test proportional hazars assumption. testing for a non-zero slope in a generalized linear regression of the scaled Schoenfeld residuals on functions of time
time.coxProp <- cox.zph(coxObj, transform = 'log')
plot(time.coxProp)
abline(h=0, lty=3, col='red')

# Make a dummy dataframe to allow Coxph lot. Take the mean age of the 2 groups
dummy <- data.frame(subtype=c('CD44','CD133'), 
                    CDE_DxAge=c(mean(boundData$CDE_DxAge[boundData$subtype %in% 'CD44'], na.rm=T), 
                           mean(boundData$CDE_DxAge[boundData$subtype %in% 'CD133'], na.rm=T)))

plot(survfit(coxObj, newdata=dummy), xlab="Survival Time (days)", ylab="Survival Probability", col=c("red",'blue'), 
main='FACS marker coexpression signature in \nGlioblastoma multiforme adjusted for age')
legend('topright', c('CD133', 'CD44'), col=c("red",'blue'), lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)
# text(locator(1),labels='p = 0.053', cex=1) #add the p-value to the graph

coxObjTmz = coxph(data.surv ~ CDE_DxAge + CDE_chemo_adjuvant_tmz + subtype , data=boundData, na.action="na.omit")
summary(coxObjTmz)

############################################# Throw in many covariates ##################################
coxObjL = coxph(data.surv ~ CDE_DxAge + subtype + G_CIMP_STATUS , data=boundData, na.action="na.omit")
summary(coxObjL)

coxObj
coxObjTmz
coxObjL

############################## Try the Agilent dataset for more n ###########################
library(GSVA)
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/clinical_dataDots.txt", row.names=1)
agilent = read.delim('~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt', row.names=1)
matched = intersect(row.names(clinical), row.names(agilent))

####### Call subtypes again ########
cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)
bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))
rm(cd133Sig, cd44Sig, cd15, aldh1, itag6, l1cam)

# Using ssGSEA heavily biases for CD44 subtype. Whereas for  GSVA the subtypes are more balanced
agilentM = t(as.matrix(agilent))
# Change NAs to 0 as GSVA doesn't like it
agilentM[is.na(agilentM)] <- 0

# Call the subtypes with GSVA
bigResult = gsva(agilentM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
bigResult = t(bigResult$es.obs)

# Merge Agilent - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
verhaakSubtype = bindGeneExprCIMPClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')

# Merge RNAseq - FACs data and clinicial data. Add Verhaak subtype
signatures = names(bigSigs)
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]
verhaakSubtype = bindGeneExprClinical(clin, bigResult, signatures)
verhaakSubtype = sort.dataframe(verhaakSubtype, 'colours')
# write.table(verhaakSubtype, '~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140909_agilentSubtypeCall.txt', sep='\t')
# Subset clinical data for intersect
# Call the FACS subtype
verhaakSubtypeCall = callMarkerSubtype(verhaakSubtype, 0, 0)

# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", 'CDE_chemo_adjuvant_tmz', 'CDE_chemo_tmz')]
#######################################

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)

############################################# Analysing the data for survival ##################################
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)
xtabs(~ CDE_DxAge + subtype, data=boundData)

coxObj = coxph(data.surv ~ CDE_DxAge + subtype, data=boundData, na.action="na.omit")
summary(coxObj)

# Test proportional hazars assumption. testing for a non-zero slope in a generalized linear regression of the scaled Schoenfeld residuals on functions of time
time.coxProp <- cox.zph(coxObj, transform = 'log')
plot(time.coxProp)
abline(h=0, lty=3, col='red')

# Make a dummy dataframe to allow Coxph lot. Take the mean age of the 2 groups
dummy <- data.frame(subtype=c('CD44','CD133'), 
                    CDE_DxAge=c(mean(boundData$CDE_DxAge[boundData$subtype %in% 'CD44'], na.rm=T), 
                                mean(boundData$CDE_DxAge[boundData$subtype %in% 'CD133'], na.rm=T)))

plot(survfit(coxObj, newdata=dummy), xlab="Survival Time (days)", ylab="Survival Probability", col=c("red",'blue'), 
     main='FACS marker coexpression signature in \nGlioblastoma multiforme adjusted for age')
legend('topright', c('CD133', 'CD44'), col=c("red",'blue'), lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)
# text(locator(1),labels='p = 0.053', cex=1) #add the p-value to the graph

coxObjTmz = coxph(data.surv ~ CDE_DxAge + CDE_chemo_adjuvant_tmz + subtype , data=boundData, na.action="na.omit")
summary(coxObjTmz)