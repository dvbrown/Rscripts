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