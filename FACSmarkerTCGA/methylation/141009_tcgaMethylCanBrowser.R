# Restart the methylation analysis using cancer browser numeric data
library(sqldf)
library(gplots)
library(RColorBrewer)
library(limma)
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