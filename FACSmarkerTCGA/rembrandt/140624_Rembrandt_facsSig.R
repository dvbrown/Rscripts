# Try and call subtype with GSVA
library(GSVA)
library(gplots)
library(RColorBrewer)
library(sqldf)

setwd('~/Documents/public-datasets/rembrandt/rembrandt_GBM/processedData/')
db <- dbConnect(SQLite(), dbname="140624_rembrandtGBM.sqlite")
dbListTables(db) 

myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

#### The signature IO ####
cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)
bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))

data = dbReadTable(db, "rembrandtSummarised")
clinical = dbReadTable(db, "rembrandtClinical", row.names='patient')

# Need to fix dis
newColnameData = gsub('_U133P2_CEL', '', colnames(data))
colnames(data) = newColnameData
row.names(data) = data$GeneSymbol

# Into the matrix for GSVA
dataM = data[,c(2:229)]
# Change NAs to 0 as GSVA doesn't like it
dataM[is.na(dataM)] <- 0

# Subset data for the matches
matched = intersect(colnames(dataM), row.names(clinical))
data.match = dataM[matched,]
clin.match = clinical[matched,]

#### Call the subtypes with GSVA ####
resultRembrandt = gsva(t(data.match), bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
resultRembrandt = t(resultRembrandt$es.obs)

#### Make heat map with my subtype ####
heatmap.2(t(data.match), cexRow=1.5, main="Enrichment of FACS marker signatures in Rembrandt Data", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          labRow=colnames(data.match), xlab="Samples", labCol=NA, offsetRow=c(1,1), margins=c(2,7.5))