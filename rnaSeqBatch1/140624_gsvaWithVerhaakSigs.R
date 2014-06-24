# Try and call subtype with GSVA
library(GSVA)
library(biomaRt)
library(gplots)
library(RColorBrewer)
source('~/Documents/Rscripts/annotateEnsembIDs.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')

setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')
verhaakSig = read.delim('~/Documents/public-datasets/TCGA/classficationSignature/131022_danFixedTCGAsignature.txt')
verhaakSigMung = list(verhaakSig[,1], verhaakSig[,2], verhaakSig[,3], verhaakSig[,4])
names(verhaakSigMung) = colnames(verhaakSig)

# Read in and sort the data as I will later remove duplicate ENSEMBL IDS
dataCorrected = read.delim('140213_normalisedCPM_facs.txt')
dataCorrected$total = rowSums(dataCorrected)
dataCorrected = sort.dataframe(dataCorrected, 7)

dataLog= read.delim('140213_normalisedLog_CPM_facs.txt')
dataLog$total = rowSums(dataLog)
dataLog = sort.dataframe(dataLog, 7)

# Extract and map gene symbol
geneNames = ensembl2officalSymbol(dataLog)
row.names(geneNames) = geneNames$SYMBOL
data.log = as.matrix(geneNames[,c(2:7)])

geneNames = ensembl2officalSymbol(dataCorrected)
row.names(geneNames) = geneNames$SYMBOL
data.corrected = as.matrix(geneNames[,c(2:7)])

result.corrected = gsva(data.corrected, verhaakSigMung)
result.corrected = t(result.corrected$es.obs)

result.log = gsva(data.log, verhaakSigMung)
result.log = t(result.log$es.obs)

myPalette <- colorRampPalette(c("blue", "white", "yellow"))(n = 1000)

# Make heat map with subtype
heatmap.2(result.corrected, cexRow=1.5, main="Enrichment of Verhaak signatures \n in RNAseq batch 1", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", Colv=NULL,
          labCol=colnames(result.corrected), xlab="RNAseq batch 1 raw cpm", labRow=row.names(result.corrected), offsetRow=c(1,1), margins=c(15,10))

# Make heat map with subtype
heatmap.2(result.log, cexRow=1.5, main="Enrichment of Verhaak signatures \n in RNAseq batch 1", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", Colv=NULL,
          labCol=colnames(result.log), xlab="RNAseq batch 1 log2 cpm", labRow=row.names(result.log), offsetRow=c(1,1), margins=c(15,10))