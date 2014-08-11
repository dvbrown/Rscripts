library(ggplot2)
library(survival)
library(coin)

source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

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
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=0, colour="red")

############################################## bind the clinical and subtyping info together #############################################

boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

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
# text(locator(1),labels='p=0.6488', cex=1) #add the p-value to the graph