library(ggplot2)
library(survival)
library(coin)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

# verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140529_verhaakSubtypeCD133_scores", row.names=1)
verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140530_liberalSignatureScores2SD.txt", row.names=1)
# The liberal signature score is more significant. Not having and intermediate case is also better

verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS")]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", "days_to_last_followup")]

# If the survival is NA, use the value for days to last followup
clin$survival = clin$CDE_survival_time
clin$survival[is.na(clin$survival)] = clin$days_to_last_followup[is.na(clin$survival)]

############################################## Segment the subtypes into CD133 CD44 #############################################

# Make an object to plot density gram
lattPlot = data.frame(as.numeric(c(verhaakSignature[,"CD133"], verhaakSignature[,"CD44"])), 
                      c(rep("CD133", nrow(verhaakSignature)), rep("CD44", nrow(verhaakSignature))))
colnames(lattPlot) = c('value', "signature")

ggplot(lattPlot, aes(value, fill = signature)) + geom_density(alpha = 0.2) +
    xlab("Signature score") + ylab("Density") + # Set axis labels
    ggtitle("Distribution of CD133 and \nCD44 signature scores") +  # Set title
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=0, colour="red") # The 0.125 is where I will call indeterminate

############################################## bind the clinical and subtyping info together #############################################
boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)
# # omit a case with no info
# empty = is.na(boundData$CDE_survival_time)
# boundData = boundData[!empty,]
############################################# Analysing the data for survival ##################################

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
test = surv_test(data.surv~boundData$subtype)#, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.0162', cex=1) #add the p-value to the graph

# Check the final table
# write.table(boundData, "./survival/survivalTables/140606_RNAseq_SurvivalBoundData.txt", sep='\t')
contingency = table(boundData$subtype, boundData$GeneExp_Subtype)
# Test the null hypothesis that the all events in the set are equally likely (from Chi-X distribution)
chisq.test(contingency)
############################################# Remove the G-CIMP cases and retest ##################################
boundDataSub = boundData[boundData$G_CIMP_STATUS == "NON G-CIMP",]

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundDataSub$CDE_survival_time, event=boundDataSub$X_EVENT)
sur.fit = survfit(data.surv~boundDataSub$subtype)

plot(sur.fit, main='FACS marker coexpression signature in Glioblastoma \nmultiforme by RNAseq (no G-CIMP)',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('CD133', 'CD44'),
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)

#test for a difference between curves
test = surv_test(data.surv~boundDataSub$subtype)
test
text(locator(1),labels='p=0.150', cex=1) #add the p-value to the graph

############################################# Survival curve for subtype ##################################
# Remove nosubtype cases
boundDataSub = boundData[!boundData$GeneExp_Subtype == "",]

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundDataSub$CDE_survival_time, event=boundDataSub$X_EVENT)

sur.fit = survfit(data.surv~GeneExp_Subtype, boundDataSub)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by RNAseq',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'green','orange', "blue"),xlim=c(0,1600), lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

legend('topright', c('Proneural', 'Neural', 'Mesenchymal',"Classical"), 
       col=c("red",'green','orange', "blue"),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundDataSub$GeneExp_Subtype)#, subset=boundDataSub$GeneExp_Subtype %in% c("Proneural", "Mesenchymal"))
test
text(locator(1),labels='p=0.6488', cex=1) #add the p-value to the graph


############################################# Plot a KM of CD133, CD44, Proneural and Mesenchymal ##################################

boundDataSub$subtype = boundDataSub$GeneExp_Subtype
boundDataSub = boundDataSub[boundDataSub$subtype %in% c('Proneural', 'Mesenchymal'),]

twoSignature = rbind(boundData, boundDataSub)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(twoSignature$CDE_survival_time, event=twoSignature$X_EVENT)

sur.fit = survfit(data.surv~subtype, twoSignature)

plot(sur.fit, main='Comparison of Verhaak and \nFACS marker signatures',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue','purple', "orange"),xlim=c(0,1600), cex=1.2, lwd=c(1.33,1.33,0.75,0.75),
     lty=c(1,1,2,2)) # hatch the verhaak lines


legend('topright', c('CD133', 'CD44', 'Mesenchymal',"Proneural"), 
       col=c("red",'blue','orange', "purple"), cex=0.9, bty='n', xjust=0.5, yjust=0.5,
       lwd=c(1.33,1.33,0.75,0.75), lty=c(1,1,2,2))

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~twoSignature$subtype, subset=twoSignature$subtype %in% c('CD133', 'CD44', 'Mesenchymal',"Proneural"))
test
text(locator(1),labels='p=0.6488', cex=1) #add the p-value to the graph

############################################# Test the Agilent dataset ##################################
# The signature derived from Agilent array
verhaakSigAgi = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140603_verhaakSubtypeAgilent.txt", row.names=1)
# Try CD133 only
verhaakSigAgi$subtype = ifelse(verhaakSigAgi[,"CD133"] >0, "CD133_high", "CD133_low")

clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/clinical_dataDots.txt", row.names=1)

verhaakSignature = verhaakSigAgi[,c("CD133","CD44")]

verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)
# verhaakSubtypeCall = verhaakSigAgi
# Extract the clinical data for the Agilent patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender", "G_CIMP_STATUS", "GeneExp_Subtype")]


# Make an object to plot density gram
lattPlot = data.frame(as.numeric(c(verhaakSignature[,"CD133"], verhaakSignature[,"CD44"])), 
                      c(rep("CD133", nrow(verhaakSignature)), rep("CD44", nrow(verhaakSignature))))
colnames(lattPlot) = c('value', "signature")

ggplot(lattPlot, aes(value, fill = signature)) + geom_density(alpha = 0.2) +
    xlab("Signature score") + ylab("Density") + # Set axis labels
    ggtitle("Distribution of CD133 and \nCD44 signature scores in Agilent") +  # Set title
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=0.05, colour="red") # The 0.125 is where I will call indeterminate


boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData)

plot(sur.fit, main='TCGA GBM patients stratified by FACS marker signature \n Agilent dataset',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),
     #xlim=c(0,1600), 
     cex=1.75, conf.int=F, lwd=1.33)

legend('topright', c('CD133', 'CD44'),
       col=c("red",'blue'),
       lwd=1.33, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

text(locator(1),labels='p=0.591', cex=1) #add the p-value to the graph
summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~subtype, boundData)#, subset=!boundData$subtype %in% "intermediate")
test

# There still appears to be a strong assosciation of the CD133 with Proneural and vice versea. Somehow no survival advantage
contingency = table(boundData$subtype, boundData$GeneExp_Subtype)
chisq.test(contingency)

############################################# Survival curve for subtype ##################################
# Remove nosubtype cases
boundDataSub = boundData[!boundData$GeneExp_Subtype == "",]

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundDataSub$CDE_survival_time, event=boundDataSub$X_EVENT)

sur.fit = survfit(data.surv~GeneExp_Subtype, boundDataSub)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by Agilent array',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'green','orange', "blue"),xlim=c(0,750), cex=1.75)

legend('topright', c('Proneural', 'Neural', 'Mesenchymal',"Classical"), 
       col=c("red",'green','orange', "blue"),lwd=1, cex=0.9, bty='n', xjust=0.5, yjust=0.5)
summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundDataSub$GeneExp_Subtype, subset=boundDataSub$GeneExp_Subtype %in% c("Proneural", "Mesenchymal"))
test


############################################# Survival curve with Agilent using the 3sd signature ##################################
verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140606_verhaakSubtypeAgilent_ssGSEA.txt", row.names=1)
verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS")]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

# Extract the clinical data for the Agilent patients
verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender")]


boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~subtype, boundData)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by Agilent',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = surv_test(data.surv~boundData$subtype)#, subset=!boundData$subtype %in% "intermediate")
test
# text(locator(1),labels='p=0.0151', cex=1) #add the p-value to the graph