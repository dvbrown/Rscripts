library(edgeR)
library(gplots)
library(RColorBrewer)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')

#Draw a heatmap. First normalise the raw counts by library size
dat = read.delim('140213_normalisedLog_CPM_facs.txt')
datM = as.matrix(dat)
myPalette <- colorRampPalette(c("white", "yellow", "red"))(n = 1000)

# Subset data for testing
datS = datM[c(1:200),]

# Extract the top 2000 most variable probes

# Make heat map with Veerhaak subtype
heatmap.2(datS, cexRow=1.5, main="", #scale='row',
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", 
          labRow=NA, xlab="Samples", labCol=colnames(datS),
          offsetRow=c(1,1), margins=c(2,7.5))