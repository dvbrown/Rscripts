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
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", "days_to_last_followup", "CDE_chemo_adjuvant_tmz", "CDE_radiation_adjuvant")]

# If the survival is NA, use the value for days to last followup
# clin$survival = clin$CDE_survival_time
# clin$survival[is.na(clin$survival)] = clin$days_to_last_followup[is.na(clin$survival)]

############################################## Segment the subtypes into CD133 CD44 #############################################

# Make an object to plot density gram
lattPlot = data.frame(as.numeric(c(verhaakSignature[,"CD133"], verhaakSignature[,"CD44"])), 
                      c(rep("CD133", nrow(verhaakSignature)), rep("CD44", nrow(verhaakSignature))))
colnames(lattPlot) = c('value', "signature")

ggplot(lattPlot, aes(value, fill = signature)) + geom_density(alpha = 0.2) +
    xlab("Signature score") + ylab("Density") + # Set axis labels
    ggtitle("Distribution of CD133 and \nCD44 signature scores") +  # Set title
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=0, colour="red")

############################################## bind the clinical and subtyping info together #############################################

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

############################################# Subset the CD133 cases and then group by TMZ exposure ##################################

cd133Patients = boundData[boundData$subtype %in% "CD133",]

############################################# Analysing the data for survival ##################################

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(cd133Patients$CDE_survival_time, event=cd133Patients$X_EVENT)

sur.fit = survfit(data.surv~CDE_chemo_adjuvant_tmz, cd133Patients)

plot(sur.fit, main='TCGA GBM cohort CD133 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('FALSE', 'TRUE'), title="Adjuvant temozolomide",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~as.factor(cd133Patients$CDE_chemo_adjuvant_tmz))#, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.137', cex=1) #add the p-value to the graph

############################################# Do the same for CD44 ##################################

cd44Patients = boundData[boundData$subtype %in% "CD44",]

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(cd44Patients$CDE_survival_time, event=cd44Patients$X_EVENT)

sur.fit = survfit(data.surv~CDE_chemo_adjuvant_tmz, cd44Patients)

plot(sur.fit, main='TCGA GBM cohort CD44 patients classified by treatment',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('FALSE', 'TRUE'), title="Adjuvant temozolomide",
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~as.factor(cd44Patients$CDE_chemo_adjuvant_tmz))#, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.0165', cex=1) #add the p-value to the graph