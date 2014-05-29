library(ggplot2)
source("~/Documents/Rscripts/120704-sortDataFrame.R")
setwd("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/")

verhaakSubtypeCall = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/survival/140529_verhaakSubtypeCD133_scores", row.names=1)
clinical = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/clinical_dataDots.txt", row.names=1)

# Extract the clinical data for the RNAseq patients
matched = intersect(row.names(clinical), row.names(verhaakSubtypeCall))
# Subset clinical data for intersect
clin = clinical[matched, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status","X_EVENT","days_to_tumor_progression", "gender")]

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
