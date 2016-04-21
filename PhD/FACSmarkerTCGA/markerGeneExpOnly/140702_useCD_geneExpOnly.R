#See what overlap just the expression CD133, CD44 and CD133 alone have with Verhaak subtype
library(sqldf)
library(gplots)
library(RColorBrewer)

getwd()
source('/Users/d.brown6/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')
setwd('~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/signatureComparison/useFACSgeneOnly/')
list.files()

#### IO ####
db <- dbConnect(SQLite(), dbname="~/Documents/public-datasets/cancerBrowser/tcgaData.sqlite")
dbListTables(db)                 # The tables in the database
dbListFields(db, "AgilentGem")       # The columns in a table

agilentData = dbReadTable(db, "AgilentGem", row.names="row_names")
clinical = dbReadTable(db, "clinicalData", row.names="row_names")
clinical = read.delim('~/Documents/public-datasets/cancerBrowser/TCGA_GBM_G4502A_07_2-2014-05-02/clinical_dataDots.txt', row.names=1)
clinPatients = gsub('\\.', '_', row.names(clinical))
row.names(clinical) = clinPatients

##### Subset the genes of interest from the gem ####
cd133 = agilentData['PROM1',]
cd44 = agilentData['CD44',]
cd15 = agilentData["FUT4",]
markerGem = t(rbind(cd133, cd44, cd15))
colnames(markerGem) = c("CD133", "CD44", "CD15")

# Bind the clincial and gene expression data
clin = clinical[, c("CDE_DxAge", "CDE_survival_time", "CDE_vital_status",
                           "G_CIMP_STATUS","GeneExp_Subtype", "X_EVENT","days_to_tumor_progression", "gender")]

data = merge(markerGem, clin, by.x="row.names", by.y='row.names')
data$colours = "black"
data$colours[data$GeneExp_Subtype == "Proneural"] = "red"
data$colours[data$GeneExp_Subtype == "Neural"] = "green"
data$colours[data$GeneExp_Subtype == "Classical"] = "blue"
data$colours[data$GeneExp_Subtype == "Mesenchymal"] = "orange"
data = sort.dataframe(data, 13)

dataM = t(as.matrix(data[,c(2,3,4)]))

myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

heatmap.2(dataM, cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", 
          Colv=data$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(data$colours), labRow=row.names(dataM), xlab="Patients", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5),
          scale='row')

heatmap.2(dataM, cexRow=1.5, main="Enrichment of FACS marker signatures \n in Molecular Subtype and G-CIMP", 
          Colv=data$colours, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(data$colours), labRow=row.names(dataM), xlab="Patients", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5),
          scale='none')
