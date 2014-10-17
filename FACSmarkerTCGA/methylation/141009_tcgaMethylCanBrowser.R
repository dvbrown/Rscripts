# Restart the methylation analysis using cancer browser numeric data
library(sqldf)
library(gplots)
library(RColorBrewer)
library(limma)
library(lumi)
source("~/Documents/Rscripts/120704-sortDataFrame.R")
source("~/Documents/Rscripts/multiplot.R")
setwd("~/Documents/public-datasets/cancerBrowser/methylation/")
list.files()

# M values says Eric

takeUnison <- function (clinicalData, genomicData) {
    # Takes the union of cases in clinical and genomic data.
    # Sort both datasets so the order is the same
  samples = intersect(colnames(genomicData), row.names(clinicalData))
  clinIntersect = clinicalData[samples,]
  gemIntersect = genomicData[,samples]
  clinIntersect <- clinIntersect[order(row.names(clinIntersect)),]
  gemIntersect <- gemIntersect[order(row.names(gemIntersect)),]
  # Returns a list with clinical data first and genomic data second
  result = list(clinIntersect, gemIntersect)
}

# Open a connection to database
db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)

###################################### Bind the clinical data and subtype calls of Agielnt and RNAseq together ##########################
# clinicalRNAseq = dbReadTable(db, "subtypeCallRNAseq")
# clinicalAgilent = dbReadTable(db, "subtypeCallAgilent")
# clinicalRNAseq = clinicalRNAseq[,c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status", "X_EVENT", "gender", "GeneExp_Subtype", "G_CIMP_STATUS", "subtype")]
# clinicalAgilent = clinicalAgilent[,c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status", "X_EVENT", "gender", "GeneExp_Subtype", "G_CIMP_STATUS", "subtype")]
# colnames(clinicalAgilent)
# colnames(clinicalRNAseq)
# clinical = rbind(clinicalAgilent, clinicalRNAseq)
# dbWriteTable(conn = db, name = "clinicalAllPatients", value = clinical, row.names = TRUE)

clinical = dbReadTable(db, "clinicalAllPatients")
# Start with h27k
h27 = read.delim("TCGA_GBM_hMethyl27-2014-08-22/genomicMatrix", row.names =1)
head(h27)

###################################### Extract samples that are in common ##########################
result = takeUnison(clinical, h27)
h27Union = result[[2]]
clinicalUnion = result[[1]]

###################################### K means clustering ##########################
# Have a look at the data distribution
methylMat = as.matrix(h27Union)
hist(methylMat)
# Want to use B value or M value or whatever lloks more normal


# Extract the 370 of the most variant probes by standard deviation
sdProbes = apply(methylMat, 1, sd)
names(sdProbes) = row.names(methylMat)
sdProbes = sort.int(sdProbes, decreasing=F)
head(sdProbes, 25)
subsetting = tail(sdProbes, 370)
methylVariable = t(methylMat[names(subsetting),])

# K-Means Cluster Analysis. Do this by patients!
fit <- kmeans(methylVariable, 3) # 3 cluster solution
# get cluster means
aggregate(methylVariable,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(methylVariable, fit$cluster)

# Merge the coexpression subtype with k means
kMclusters = merge.data.frame(mydata[,c(370,371)], clinicalUnion, by.x="row.names", by.y="row.names")
    
###################################### limma analysis ##########################
f = factor(clinicalUnion$subtype)
design = model.matrix(~f)
colnames(design) = levels(f)
fit = lmFit(h27Union, design)

cont.matrix = makeContrasts(cd133vscd44="CD44", levels=design)
fit2  = contrasts.fit(fit, cont.matrix)
fit2  = eBayes(fit2)

#write the output to a table of differentially expressed genes. Change this value to suit
result = topTable(fit2, coef=1, number=22277, #genelist=genelist,  
                  adjust='BH', sort.by='p', lfc=0)

head(result,100)
results <- decideTests(fit2)
summary(results)
result$threshold = as.factor(abs(result$logFC) > 2 & result$adj.P.Val < 0.1)

###################################### Remove G-CIMP cases and retest ##########################
clinicalnoCIMP = clinicalUnion[clinicalUnion$G_CIMP_STATUS %in% "NON G-CIMP",]
result2 = takeUnison(clinicalnoCIMP, h27)
h27Union = result2[[2]]
clinicalUnion = result2[[1]]

f = factor(clinicalUnion$subtype)
design = model.matrix(~f)
colnames(design) = levels(f)
fit = lmFit(h27Union, design)

cont.matrix = makeContrasts(cd133vscd44="CD44", levels=design)
fit2  = contrasts.fit(fit, cont.matrix)
fit2  = eBayes(fit2)

#write the output to a table of differentially expressed genes. Change this value to suit
result2 = topTable(fit2, coef=1, number=22277, #genelist=genelist,  
                  adjust='BH', sort.by='logFC', lfc=0)

head(result2,100)
result3 <- decideTests(fit2)
summary(result3)

##Construct the volcano plot object.
##ighlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
result2$threshold = as.factor(abs(result2$logFC) > 2 & result2$adj.P.Val < 0.1)

volcanoNoCIMP = ggplot(data=result2, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) + scale_colour_manual(values=c('grey', 'blue')) +
    xlim(c(-0.2, 0.2)) + ylim(c(0, 0.67)) +
    xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value") + # Set axis labels
    ggtitle("TCGA GBM dataset Infinium HumanMethylation27\n G-CIMP cases removed") +  # Set title
    theme_bw(base_size=18)

volcanoCIMP = ggplot(data=result, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) + scale_colour_manual(values=c('grey', 'blue')) +
    xlim(c(-0.2, 0.2)) + ylim(c(0, 0.67)) +
    xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value") + # Set axis labels
    ggtitle("TCGA GBM dataset Infinium HumanMethylation27\n G-CIMP cases retained") +  # Set title
    theme_bw(base_size=18)

# In conclusion there is a difference but the magnitude is not large < 0.5 lfc
# If you remove the G-CIMP cases then there is no significance

multiplot(volcanoCIMP, volcanoNoCIMP, cols=2)