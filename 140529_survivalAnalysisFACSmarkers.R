library(ggplot2)
library(survival)

source("~/Documents/Rscripts/120704-sortDataFrame.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

callMarkerSubtype <- function (signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # Not having and intermediate case is also better for the Kaplan Myer curve
    #signatureScore$subtype[signatureScore[,"CD133"] < -CD133cutoff & signatureScore[,"CD44"] < CD44cutoff] = "intermediate"
    signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}

# verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140529_verhaakSubtypeCD133_scores", row.names=1)
verhaakSignature = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140530_liberalSignatureScores2SD.txt", row.names=1)
# The liberal signature score is more significant. Not having and intermediate case is also better

verhaakSignature = verhaakSignature[,c("CD133","CD44","GeneExp_Subtype","G_CIMP_STATUS")]
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender")]

############################################## Segment the subtypes into CD133 CD44 and indetermiant #############################################

# Make an object to plot density gram
lattPlot = data.frame(as.numeric(c(verhaakSignature[,"CD133"], verhaakSignature[,"CD44"])), 
                      c(rep("CD133", nrow(verhaakSignature)), rep("CD44", nrow(verhaakSignature))))
colnames(lattPlot) = c('value', "signature")

ggplot(lattPlot, aes(value, fill = signature)) + geom_density(alpha = 0.2) +
    xlab("Signature score") + ylab("Density") + # Set axis labels
    ggtitle("Distribution of CD133 and \nCD44 signature scores") +  # Set title
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=0, colour="red") # The 0.125 is where I will call indeterminate

verhaakSubtypeCall = callMarkerSubtype(verhaakSignature, 0, 0)

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

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by RNAseq',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~boundData$subtype)#, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.0151', cex=1) #add the p-value to the graph




############################################# Remove the G-CIMP cases and retest ##################################
boundDataSub = boundData[boundData$G_CIMP_STATUS == "NON G-CIMP",]

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundDataSub$CDE_survival_time, event=boundDataSub$X_EVENT)
sur.fit = survfit(data.surv~boundDataSub$subtype)

plot(sur.fit, main='FACS marker coexpression signature in Glioblastoma \nmultiforme by RNAseq no (G-CIMP)',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue'),#'green'),
     xlim=c(0,750), cex=1.75, conf.int=F, lwd=1.5)

legend('topright', c('CD133', 'CD44'),# 'Intermediate'), 
       col=c("red",'blue'),#'green'),
       lwd=2, cex=1.2, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)

#test for a difference between curves
test = survdiff(data.surv~boundDataSub$subtype, subset=!boundDataSub$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.152', cex=1) #add the p-value to the graph

############################################# Survival curve for subtype ##################################
#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)

sur.fit = survfit(data.surv~GeneExp_Subtype, boundData)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by RNAseq',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'green','orange', "blue"),xlim=c(0,750), cex=1.75)

legend('topright', c('Proneural', 'Neural', 'Mesenchymal',"Classical"), 
       col=c("red",'green','orange', "blue"),lwd=1, cex=0.9, bty='n', xjust=0.5, yjust=0.5)
summary(data.surv)
#test for a difference between curves
test = survdiff(data.surv~boundData$GeneExp_Subtype, subset=boundData$GeneExp_Subtype %in% c("Proneural", "Mesenchymal"))
test
# text(locator(1),labels='p=0.071', cex=1) #add the p-value to the graph