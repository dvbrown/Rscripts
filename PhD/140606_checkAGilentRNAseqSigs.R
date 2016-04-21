# Compare the agilent array and the RNAseq data
source("~/Documents/Rscripts/140211_multiplotGgplot2.R")

rnaseq = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140530_liberalSignatureScores2SD.txt", row.names=1)
aglient = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140603_verhaakSubtypeAgilent.txt", row.names=1)
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

rnaseq$origin = "RNA"
rnaseq = rnaseq[,c("CD133", "CD44", "GeneExp_Subtype", "G_CIMP_STATUS", "origin")]
aglient$origin = "Array"
aglient = aglient[,c("CD133", "CD44", "GeneExp_Subtype", "G_CIMP_STATUS", "origin")]

lattPlot = rbind(rnaseq, aglient)

a = ggplot(lattPlot, aes(CD133, fill = origin)) + geom_density(alpha = 0.2) +
    xlab("Signature score") + ylab("Density") + # Set axis labels
    ggtitle("Distribution of CD133 and signature scores in patient data") +  # Set title
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=0.0, colour="red") # The 0.125 is where I will call indeterminate

b = ggplot(lattPlot, aes(CD44, fill = origin)) + geom_density(alpha = 0.2) +
    xlab("Signature score") + ylab("Density") + # Set axis labels
    ggtitle("Distribution of CD44 signature scores in patient data") +  # Set title
    coord_cartesian(xlim = c(-1, 1)) + theme_bw(base_size=20) + geom_vline(xintercept=0.0, colour="red") # The 0.125 is where I will call indeterminate

multiplot(a,b)

lattPlot$subtype = ifelse(lattPlot[,"CD133"] > lattPlot[,"CD44"], "CD133", "CD44")
lattPlot$subtype[lattPlot[,"CD133"] < 0  & lattPlot[,"CD44"] < 0] = "intermediate"
lattPlot$subtype = as.factor(lattPlot$subtype)

contingency = table(lattPlot$origin, lattPlot$subtype)
fisher.test(contingency)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(lattPlot))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT", "gender")]

boundData = merge.data.frame(clin, lattPlot, by.x="row.names", by.y="row.names")
boundData = sort.dataframe(boundData, "subtype")
row.names(boundData) = boundData$Row.names
boundData$subtype = as.factor(boundData$subtype)
boundData$gender = as.factor(boundData$gender)

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
test = surv_test(data.surv~boundData$subtype)#, subset=!boundData$subtype %in% "intermediate")
test
text(locator(1),labels='p=0.0151', cex=1) #add the p-value to the graph
