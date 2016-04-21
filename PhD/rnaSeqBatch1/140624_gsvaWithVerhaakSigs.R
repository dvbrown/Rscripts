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
result.corrected = (result.corrected$es.obs)

result.log = gsva(data.log, verhaakSigMung)
result.log = result.log$es.obs
colnames(result.corrected) = c("MU011", "MU020", "MU034", "MU035", "MU039", "MU041")

myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

# Make heat map with subtype
heatmap.2(result.corrected, cexRow=1.5, main="", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", #Colv=NULL,
          labCol=colnames(result.corrected), xlab="PDGC sample", labRow=row.names(result.corrected), 
          Rowv =c("Proneural", "Neural", "Classical", "Mesenchymal"), offsetRow=c(1,1), margins=c(8, 10))

# Make heat map with subtype
# heatmap.2(result.log, cexRow=1.5, main="Enrichment of Verhaak signatures \n in RNAseq batch 1", 
#           keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", Colv=NULL,
#           labCol=colnames(result.log), xlab="RNAseq batch 1 log2 cpm", labRow=row.names(result.log), offsetRow=c(1,1), margins=c(15,10))


################################# Enrich against my CD133/ CD44 signature ######################################
cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15))

result.log = gsva(data.log, bigSigs, verbose=T, parallel.sz=1)
result.log = (result.log$es.obs)

# Make heat map with subtype
heatmap.2(result.log, cexRow=1.5, main="Enrichment of CD133/ CD44 signatures \n in RNAseq batch 1", 
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", Colv=NULL,
          labCol=colnames(result.log), xlab="RNAseq batch 1 log2 cpm", labRow=row.names(result.log), offsetRow=c(1,1), margins=c(15,10))


biggerSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
                  'Mesenchymal'=verhaakSigMung$Mesenchymal, 'Proneural'= verhaakSigMung$Proneural)

resultBig.log = gsva(data.log, biggerSigs, verbose=T, parallel.sz=1)
resultBig.log = resultBig.log$es.obs

# Make heat map with subtype
# heatmap.2(resultBig.log, cexRow=1.5, main="Enrichment of signatures \n in RNAseq batch 1", 
#           keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", Colv=NULL,
#           labCol=colnames(resultBig.log), xlab="RNAseq batch 1 log2 cpm", labRow=row.names(resultBig.log), offsetRow=c(1,1), margins=c(15,10))