library(ggplot2)
library(survival)
library(coin)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140603_verhaakSubtypeAgilent.txt", row.names=1)

verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS")]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", 'CDE_chemo_adjuvant_tmz', 'CDE_chemo_tmz',
                           'CDE_radiation_any', 'CDE_tmz_chemoradiation_standard', 'GeneExp_Subtype', 'G_CIMP_STATUS')]

# Remove cases where NA for treatment
#clin = clin[!is.na(clin$CDE_chemo_tmz),]

############################################## bind the clinical and subtyping info together #############################################
boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

# Fix up the GCIMP to only have true and false
boundData$G_CIMP_STATUS = boundData$G_CIMP_STATUS.x
boundData$G_CIMP_STATUS = as.character(boundData$G_CIMP_STATUS)
boundData$G_CIMP_STATUS[boundData$G_CIMP_STATUS %in% 'G-CIMP'] = TRUE
boundData$G_CIMP_STATUS[!boundData$G_CIMP_STATUS %in% TRUE] = FALSE
boundData$G_CIMP_STATUS = as.character(boundData$G_CIMP_STATUS)
boundData$G_CIMP_STATUS = as.logical(boundData$G_CIMP_STATUS)
boundData$CDE_radiation_any = as.logical(boundData$CDE_radiation_any)
boundData = boundData[,c(1:11,13,14,17,18)]

############################################# Subset the Subtype cases indivdually ##################################

cd133Patients = boundData[boundData$subtype %in% "CD133",]
cd44Patients = boundData[boundData$subtype %in% "CD44",]

############################################# Analysing the data for survival ##################################

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

#generate the survival object and plot a Kaplan-Meier
data.surv.cd133 = Surv(cd133Patients$CDE_survival_time, event=cd133Patients$X_EVENT)
data.surv.cd44 = Surv(cd44Patients$CDE_survival_time, event=cd44Patients$X_EVENT)

############################# Make a Cox proportional hazards model to get hazard ratios ##############################
# CD133 patients
coxCD133ph = coxph(data.surv.cd133 ~ CDE_DxAge + G_CIMP_STATUS + CDE_chemo_tmz + CDE_radiation_any, 
                   data=cd133Patients, na.action="na.omit")

coxCD44ph = coxph(data.surv.cd44 ~ CDE_DxAge + G_CIMP_STATUS + CDE_chemo_tmz + CDE_radiation_any, 
                  data=cd44Patients, na.action="na.omit")
summary(coxCD133ph)
summary(coxCD44ph)




# par(mfrow=c(2,2))
par(mfrow=c(1,1))
############################# CD133 Temozolomide ##############################
attach(cd133Patients)
cd133.temo = data.frame(CDE_chemo_tmz=c(TRUE, FALSE), CDE_DxAge=rep(mean(CDE_DxAge, na.rm=T),2), G_CIMP_STATUS=rep(1-mean(G_CIMP_STATUS),2),
                        CDE_radiation_any=rep(mean(CDE_radiation_any, na.rm=T),2))
detach(cd133Patients)

plot(survfit(coxCD133ph, newdata=cd133.temo, na.action=na.pass), main='TCGA GBM cohort CD133 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'), cex=1.75, conf.int=F, lwd=1.33)
legend('topright', c('FALSE', 'TRUE'), title="Temozolomide",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

############################# CD133 Radiation ##############################
attach(cd133Patients)
cd133.rad = data.frame(CDE_radiation_any=c(TRUE, FALSE), CDE_DxAge=rep(mean(CDE_DxAge, na.rm=T),2), G_CIMP_STATUS=rep(1-mean(G_CIMP_STATUS),2),
                       CDE_chemo_tmz=rep(mean(CDE_chemo_tmz, na.rm=T),2))
detach(cd133Patients)

plot(survfit(coxCD133ph, newdata=cd133.rad, na.action=na.pass), main='TCGA GBM cohort CD133 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'), cex=1.75, conf.int=F, lwd=1.33)
legend('topright', c('FALSE', 'TRUE'), title="Radiation",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)


############################# CD44 Temozolomide ##############################
attach(cd44Patients)
cd44.temo = data.frame(CDE_chemo_tmz=c(TRUE, FALSE), CDE_DxAge=rep(mean(CDE_DxAge, na.rm=T),2), G_CIMP_STATUS=rep(1-mean(G_CIMP_STATUS),2),
                        CDE_radiation_any=rep(mean(CDE_radiation_any, na.rm=T),2))
detach(cd44Patients)

plot(survfit(coxCD44ph, newdata=cd44.temo, na.action=na.pass), main='TCGA GBM cohort CD44 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'), cex=1.75, conf.int=F, lwd=1.33)
legend('topright', c('FALSE', 'TRUE'), title="Temozolomide",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

############################# CD44 Radiation ##############################
attach(cd44Patients)
cd44.rad = data.frame(CDE_radiation_any=c(TRUE, FALSE), CDE_DxAge=rep(mean(CDE_DxAge, na.rm=T),2), G_CIMP_STATUS=rep(1-mean(G_CIMP_STATUS),2),
                       CDE_chemo_tmz=rep(mean(CDE_chemo_tmz, na.rm=T),2))
detach(cd44Patients)

plot(survfit(coxCD44ph, newdata=cd44.rad, na.action=na.pass), main='TCGA GBM cohort CD44 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'), cex=1.75, conf.int=F, lwd=1.33)
legend('topright', c('FALSE', 'TRUE'), title="Radiation",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)