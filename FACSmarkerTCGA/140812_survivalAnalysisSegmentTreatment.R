library(ggplot2)
library(survival)
library(coin)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140603_verhaakSubtypeAgilent.txt", row.names=1)
# Use the Agilent dataset as this has more cases 

verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS")]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", 'CDE_chemo_adjuvant_tmz', 'CDE_chemo_tmz',
                           'CDE_radiation_any', 'CDE_tmz_chemoradiation_standard', 'GeneExp_Subtype', 'G_CIMP_STATUS')]

############################################## bind the clinical and subtyping info together #############################################

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

# Fix up the GCIMP to only have true and false
boundData$G_CIMP_STATUS = as.character(boundData$G_CIMP_STATUS)
boundData$G_CIMP_STATUS[boundData$G_CIMP_STATUS %in% 'G-CIMP'] = TRUE
boundData$G_CIMP_STATUS[!boundData$G_CIMP_STATUS %in% 'TRUE'] = FALSE
boundData$G_CIMP_STATUS = as.factor(boundData$G_CIMP_STATUS)

############################################# Subset the Subtype cases indivdually ##################################

cd133Patients = boundData[boundData$subtype %in% "CD133",]
cd44Patients = boundData[boundData$subtype %in% "CD44",]

############################################# Analysing the data for survival ##################################

#generate the survival object and plot a Kaplan-Meier
data.surv.cd133 = Surv(cd133Patients$CDE_survival_time, event=cd133Patients$X_EVENT)
data.surv.cd44 = Surv(cd44Patients$CDE_survival_time, event=cd44Patients$X_EVENT)

# sur.fit.cd133 = survfit(data.surv.cd133~CDE_radiation_any, cd133Patients)
# sur.fit.cd44 = survfit(data.surv.cd44~CDE_radiation_any, cd44Patients)
sur.fit.cd133 = survfit(data.surv.cd133~CDE_chemo_tmz, cd133Patients)
sur.fit.cd44 = survfit(data.surv.cd44~CDE_chemo_tmz, cd44Patients)

par(mfrow=c(2,1))
plot(sur.fit.cd133, main='TCGA GBM cohort CD133 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     cex=1.75, conf.int=F, lwd=1.33)
legend('topright', c('FALSE', 'TRUE'), title="Adjuvant temozolomide",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv.cd133)
test.cd133 = surv_test(data.surv.cd133~as.factor(cd133Patients$CDE_chemo_tmz))
test.cd133
# text(locator(1),labels='p=0.137', cex=1) #add the p-value to the graph

plot(sur.fit.cd44, main='TCGA GBM cohort CD44 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     cex=1.75, conf.int=F, lwd=1.33)
legend('topright', c('FALSE', 'TRUE'), title="Adjuvant temozolomide",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)
summary(data.surv.cd44)
#test for a difference between curves
test.cd44 = surv_test(data.surv.cd44~as.factor(cd44Patients$CDE_chemo_tmz))
test.cd44
# text(locator(1),labels='p=0.0165', cex=1) #add the p-value to the graph

data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)
sur.fit = survfit(data.surv ~ subtype, boundData)
par(mfrow=c(1,1))
plot(sur.fit, main='TCGA GBM cohort all patients classified by subtype',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     cex=1.75, conf.int=F, lwd=1.33)
legend('topright', c('CD133', 'CD44'), title="Coexpression subtype",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
test = surv_test(data.surv~as.factor(boundData$subtype))
test
##########################################################################################

#### Look at radiation ####
sur.fit = survfit(data.surv~CDE_radiation_adjuvant, cd133Patients)

plot(sur.fit, main='TCGA GBM cohort CD133 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('FALSE', 'TRUE'), title="Adjuvant radiation",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~as.factor(cd133Patients$CDE_radiation_adjuvant))#, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=1.4e-06', cex=1) #add the p-value to the graph



#### CD44 and radiation ####
data.surv = Surv(cd44Patients$CDE_survival_time, event=cd44Patients$X_EVENT)
sur.fit = survfit(data.surv~CDE_radiation_adjuvant, cd44Patients)

plot(sur.fit, main='TCGA GBM cohort CD44 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('FALSE', 'TRUE'), title="Adjuvant radiation",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~as.factor(cd44Patients$CDE_radiation_adjuvant))#, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.004', cex=1) #add the p-value to the graph