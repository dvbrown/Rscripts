library(ggplot2)
library(survival)
library(coin)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

############################################## IO and Munging #############################################

verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140530_liberalSignatureScores2SD.txt", row.names=1)

verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS")]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", "days_to_last_followup")]

# If the survival is NA, use the value for days to last followup. THIS IS DIFFERENT
clin$survival = clin$CDE_survival_time
clin$survival[is.na(clin$survival)] = clin$days_to_last_followup[is.na(clin$survival)]

############################################## bind the clinical and subtyping info together #####################################

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

############################################# Analysing the data for survival the original way ##################################

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData)

plot(sur.fit, main='TCGA GBM cohort classified by FACS marker signature',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('CD133', 'CD44'),
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype)
test
# text(locator(1),labels='p=0.0162', cex=1) #add the p-value to the graph

############################################# Analysing the data for survival the correct consoring way ##################################

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$survival, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData)

plot(sur.fit, main='TCGA GBM cohort classified by FACS marker signature',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('CD133', 'CD44'),
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype)
test
# text(locator(1),labels='p=0.0696', cex=1) #add the p-value to the graph

############################################# Analysing the data for survival the correct consoring way. Use vital status this time ##################################
# THIS IS EQUIVALENT TO THE ORIGINAL WAY

# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", "days_to_last_followup")]

# If the survival is NA, use the value for days to last followup. THIS IS DIFFERENT
clin$survival = clin$CDE_survival_time
clin$survival[is.na(clin$survival)] = clin$days_to_last_followup[is.na(clin$survival)]

# Use vital status as event
clin$status = NA
clin$status[clin$CDE_vital_status %in% "DECEASED"] = 1
clin$status[clin$CDE_vital_status %in% "LIVING"] = 0

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$survival, event=boundData$status)

sur.fit = survfit(data.surv~subtype, boundData)

plot(sur.fit, main='TCGA GBM cohort classified by FACS marker signature',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('CD133', 'CD44'),
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype)
test