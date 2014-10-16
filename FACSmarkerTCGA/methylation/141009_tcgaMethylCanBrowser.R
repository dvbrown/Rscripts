# Restart the methylation analysis using cancer browser numeric data
library(sqldf)
library(gplots)
library(RColorBrewer)
library(limma)
source("~/Documents/Rscripts/120704-sortDataFrame.R")
setwd("~/Documents/public-datasets/cancerBrowser/methylation/")
list.files()

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

f = factor(clinicalUnion$subtype)
design = model.matrix(~f)
colnames(design) = levels(f)
fit = lmFit(h27Union, design)

cont.matrix = makeContrasts(cd133vscd44="CD44", levels=design)
fit2  = contrasts.fit(fit, cont.matrix)
fit2  = eBayes(fit2)

#write the output to a table of differentially expressed genes. Change this value to suit
result = topTable(fit2, coef=1, number=22277, #genelist=genelist,  
                  adjust='BH', sort.by='logFC', lfc=0)

head(result,100)
results <- decideTests(fit2)
summary(results)
result$threshold = as.factor(abs(result$logFC) > 2 & result$adj.P.Val < 0.05)

volcanoCIMP = ggplot(data=result, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) + scale_colour_manual(values='grey') +
    xlim(c(-0.2, 0.2)) + ylim(c(0, 0.3)) +
    xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value") + # Set axis labels
    ggtitle("TCGA GBM dataset Infinium HumanMethylation27\n G-CIMP cases retained") +  # Set title
    theme_bw(base_size=18)
volcanoCIMP

###################################### Remove G-CIMP cases and retest ##########################
clinicalnoCIMP = clinicalUnion[clinicalUnion$G_CIMP_STATUS %in% "NON G-CIMP",]
result = takeUnison(clinicalnoCIMP, h27)
h27Union = result[[2]]
clinicalUnion = result[[1]]

f = factor(clinicalUnion$subtype)
design = model.matrix(~f)
colnames(design) = levels(f)
fit = lmFit(h27Union, design)

cont.matrix = makeContrasts(cd133vscd44="CD44", levels=design)
fit2  = contrasts.fit(fit, cont.matrix)
fit2  = eBayes(fit2)

#write the output to a table of differentially expressed genes. Change this value to suit
result = topTable(fit2, coef=1, number=22277, #genelist=genelist,  
                  adjust='BH', sort.by='logFC', lfc=0)

head(result,100)
results <- decideTests(fit2)
summary(results)

##Construct the volcano plot object.
##ighlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
result$threshold = as.factor(abs(result$logFC) > 2 & result$adj.P.Val < 0.05)

volcanoNoCIMP = ggplot(data=result, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) + scale_colour_manual(values='grey') +
    xlim(c(-0.2, 0.2)) + ylim(c(0, 0.3)) +
    xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value") + # Set axis labels
    ggtitle("TCGA GBM dataset Infinium HumanMethylation27\n G-CIMP cases removed") +  # Set title
    theme_bw(base_size=18)
volcanoNoCIMP

# In conclusion there is a difference but the magnitude is not large < 0.5 lfc
# If you remove the G-CIMP cases then there is no significance