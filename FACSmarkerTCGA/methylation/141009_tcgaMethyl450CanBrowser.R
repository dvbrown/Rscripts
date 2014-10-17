# Analyse the 450k methylation this time
library(sqldf)
library(gplots)
library(RColorBrewer)
library(limma)
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

clinical = dbReadTable(db, "clinicalAllPatients")
# Start with h27k
h450 = read.delim("TCGA_GBM_hMethyl450-2014-08-22/genomicMatrix", row.names =1)
head(h450)

###################################### Extract samples that are in common ##########################
result = takeUnison(clinical, h450)
h450Union = result[[2]]
clinicalUnion = result[[1]]

f = factor(clinicalUnion$subtype)
design = model.matrix(~f)
colnames(design) = levels(f)
fit = lmFit(h450Union, design)

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