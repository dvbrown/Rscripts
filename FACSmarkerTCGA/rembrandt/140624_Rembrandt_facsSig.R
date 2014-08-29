# Try and call subtype with GSVA
library(GSVA)
library(gplots)
library(RColorBrewer)
library(sqldf)
source('~/Documents/Rscripts/120704-sortDataFrame.R')
source("~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R")

setwd('~/Documents/public-datasets/rembrandt/rembrandt_GBM/processedData/')
db <- dbConnect(SQLite(), dbname="140624_rembrandtGBM.sqlite")
dbListTables(db) 

myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

#### The signature IO ####
cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))
rm(cd133Sig, cd44Sig, cd15, aldh1, itag6, l1cam)


verhaakSig = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')
verSigs = list(verhaakSig$Proneural, verhaakSig$Neural, verhaakSig$Classical, verhaakSig$Mesenchymal)
names(verSigs) = colnames(verhaakSig)

# The highest probe
# data = dbReadTable(db, "rembrandtSummarised", row.names="GeneSymbol")

# The mean of probes
data = dbReadTable(db, "rembrandtProbeMean", row.names="GeneSymbol")
clinical = dbReadTable(db, "rembrandtClinical", row.names='patient')

# Need to fix dis
newColnameData = gsub('_U133P2_CEL', '', colnames(data))
colnames(data) = newColnameData

hist(as.numeric(data['PROM1',]), breaks='Scott')
hist(as.numeric(data['CD44',]), breaks='Scott')

# Into the matrix for GSVA
dataM = as.matrix(data[,c(1:229)])
# Change NAs to 0 as GSVA doesn't like it
dataM[is.na(dataM)] <- 0

# Subset data for the matches
matched = intersect(colnames(dataM), row.names(clinical))
data.match = dataM[,matched]
clin.match = clinical[matched,]

#### Call the subtypes with GSVA ####
resultVerhaak = gsva(data.match, verSigs,  rnaseq=F, verbose=T, parallel.sz=1)
resultVerhaak = t(resultVerhaak$es.obs)

# Write some code to call the verhaak subtype
resultVerhaakIndex = as.data.frame(resultVerhaak)
index = max.col(resultVerhaakIndex)
resultVerhaakIndex = cbind(resultVerhaakIndex, index)
resultVerhaakIndex$subtype = ""
resultVerhaakIndex$subtype[resultVerhaakIndex$index == 1] = 'red'
resultVerhaakIndex$subtype[resultVerhaakIndex$index == 2] = 'green'
resultVerhaakIndex$subtype[resultVerhaakIndex$index == 3] = 'blue'
resultVerhaakIndex$subtype[resultVerhaakIndex$index == 4] = 'orange'
resultVerhaakIndex = sort.dataframe(resultVerhaakIndex, 5, highFirst=F)
resultVerhaak = as.matrix(resultVerhaakIndex[,c(1:4)])

resultRembrandt = gsva(data.match, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
resultRembrandt = t(resultRembrandt$es.obs)

# Get the colours of the subtype into the same order for the rembrandt called subtypes IMPORTANT!
resultRembrandtMerge = merge(resultRembrandt, resultVerhaakIndex[,c(5,6)], by.x='row.names', by.y='row.names')
row.names(resultRembrandtMerge) = resultRembrandtMerge$Row.names
resultRembrandtMerge = sort.dataframe(resultRembrandtMerge, 8, highFirst=F)
resultRembrandt = as.matrix(resultRembrandtMerge[,c(2:7)])

#### Make heat map with my subtype ####
heatmap.2(t(resultVerhaak), cexRow=1.5, main="Identifying molecular subtype in Rembrandt Data", 
          Colv=resultVerhaakIndex$subtype, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(resultVerhaakIndex$subtype), labRow=colnames(resultVerhaakIndex), xlab="Rembrandt samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

heatmap.2(t(resultRembrandt), cexRow=1.5, main="Enrichment of FACS marker signatures\nin Rembrandt GBM", 
          Colv=resultRembrandtMerge$subtype, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(resultRembrandtMerge$subtype), labRow=colnames(resultRembrandt), xlab="Rembrandt samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

# Reduced subtype
# Make heat map with only 3 markers
heatmap.2(t(resultRembrandt[,c(1:3)]), cexRow=1.5, main="Enrichment of FACS marker signatures\nin Rembrandt GBM", 
          Colv=resultRembrandtMerge$subtype, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(resultRembrandtMerge$subtype), labRow=colnames(resultRembrandt), xlab="Rembrandt samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

# Make heat map with only 3 mRNAs
markers = as.matrix(data.match[c('CD44', 'FUT4', 'PROM1'),])
markers = apply(markers, c(1,2), as.numeric)

heatmap.2(markers, cexRow=1.5, main="Enrichment of FACS marker mRNAs\n in Molecular Subtype and G-CIMP", scale='row',
          Colv=resultRembrandtMerge$subtype, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(resultRembrandtMerge$subtype), labRow=colnames(resultRembrandt), xlab="Rembrandt samples", labCol=NA, 
          offsetRow=c(1,1), margins=c(2,7.5), ylab="Marker")

#################### Call the molecualr subtype ####################

subtypeRembrandt = callMarkerSubtype(as.data.frame(resultRembrandtMerge[,c(2,3,8,9)]), 0, 0)
# subtypeRembrandt = merge(subtypeRembrandt, resultVerhaakIndex, by.x='row.names', by.y='row.names')

# dbWriteTable(conn = db, name = "facsSubtyeRembrandtProbeMean", value = resultRembrandtMerge, row.names = TRUE)
# dbWriteTable(conn = db, name = "facsSubtyeRembrandtProbeHighest", value = subtypeRembrandt, row.names = TRUE)

# Build a contingency table and test membership
subtypeRembrandt$verhaak = ""
subtypeRembrandt$verhaak[subtypeRembrandt$index == 1] = 'Proneural'
subtypeRembrandt$verhaak[subtypeRembrandt$index == 2] = 'Neural'
subtypeRembrandt$verhaak[subtypeRembrandt$index == 3] = 'Classical'
subtypeRembrandt$verhaak[subtypeRembrandt$index == 4] = 'Mesenchymal'

subtypeRembrandtPM = subtypeRembrandt[subtypeRembrandt$verhaak %in% c('Proneural', 'Mesenchymal'),]
tab = table(subtypeRembrandtPM$subtype, subtypeRembrandtPM$verhaak)
fisher.test(tab)

subtypeRembrandtCN = subtypeRembrandt[subtypeRembrandt$verhaak %in% c('Neural', 'Classical'),]
tab = table(subtypeRembrandtCN$subtype, subtypeRembrandtCN$verhaak)
fisher.test(tab)