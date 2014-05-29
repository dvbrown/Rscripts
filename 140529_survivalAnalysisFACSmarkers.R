library(ggplot2)
library(survival)

source("~/Documents/Rscripts/120704-sortDataFrame.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

buildClassifier = function(signatureSurvivalFrame, percentileNum) {
    # Take a dataframe containing a signature score and censorship status and add the group membership eg 'high' or 'low'
    # Allows one to vary the percentile used as the classifier
    percent = quantile(signatureSurvivalFrame$sigScore, probs=percentileNum, names=T)
    signatureSurvivalFrame$percentile = ifelse(signatureSurvivalFrame$sigScore >= percent, 'high', 'low')
    return (signatureSurvivalFrame)
}

verhaakSubtypeCall = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140529_verhaakSubtypeCD133_scores", row.names=1)
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender")]

############################################## Segment the subtypes into CD133 CD44 and indetermiant #############################################

# Make an object to plot density gram
lattPlot = data.frame(as.numeric(c(verhaakSubtypeCall[,"CD133"], verhaakSubtypeCall[,"CD44"])), 
                      c(rep("CD133", nrow(verhaakSubtypeCall)), rep("CD44", nrow(verhaakSubtypeCall))))
colnames(lattPlot) = c('value', "signature")

ggplot(lattPlot, aes(value, fill = signature)) + geom_density(alpha = 0.2) +
    xlab("Signature score") + ylab("Density") + # Set axis labels
    ggtitle("Distribution of CD133 and \nCD44 signature scores") +  # Set title
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=-0.125, colour="red") # The 0.125 is where I will call indeterminate

# Now call the subtypes
verhaakSubtypeCall$subtype = ""
verhaakSubtypeCall$subtype = ifelse(verhaakSubtypeCall[,"CD133"] > verhaakSubtypeCall[,"CD44"], "CD133", "CD44")
verhaakSubtypeCall$subtype[verhaakSubtypeCall$CD133 < -0.125 & verhaakSubtypeCall$CD44 < 0] = "intermediate"
verhaakSubtypeCall = sort.dataframe(verhaakSubtypeCall, "subtype")

############################################## bind the clinical and subtyping info together #############################################
boundData = merge.data.frame(clin, verhaakSubtypeCall, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
boundData$subtype = as.factor(boundData$subtype)
############################################# Analysing the data for survival ##################################

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundData$CDE_survival_time, event=boundData$X_EVENT)
sur.fit = survfit(data.surv~boundData$subtype)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by RNAseq',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue','green'),xlim=c(0,750), cex=1.75)

legend('topright', c('CD133', 'CD44', 'Intermediate'), 
       col=c("red",'blue','green'),lwd=1, cex=0.9, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)

#test for a difference between curves
test = survdiff(data.surv~boundData$subtype, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.071', cex=1) #add the p-value to the graph

############################################# Remove the G-CIMP cases and retest ##################################
boundDataSub = boundData[boundData$G_CIMP_STATUS == "NON G-CIMP",]

#generate the survival object and plot a Kaplan-Meier
data.surv = Surv(boundDataSub$CDE_survival_time, event=boundDataSub$X_EVENT)
sur.fit = survfit(data.surv~boundDataSub$subtype)

plot(sur.fit, main='FACS marker coexpression signature in \nGlioblastoma multiforme by RNAseq (no G-CIMP)',ylab='Survival probability',xlab='survival (days)', 
     col=c("red",'blue','green'),xlim=c(0,750), cex=1.75)

legend('topright', c('CD133', 'CD44', 'Intermediate'), 
       col=c("red",'blue','green'),lwd=1, cex=0.9, bty='n', xjust=0.5, yjust=0.5)

summary(data.surv)

#test for a difference between curves
test = survdiff(data.surv~boundDataSub$subtype, subset=!boundDataSub$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.31', cex=1) #add the p-value to the graph