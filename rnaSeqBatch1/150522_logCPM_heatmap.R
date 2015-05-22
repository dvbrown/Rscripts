library(edgeR)
library(gplots)
library(RColorBrewer)
source('~/Documents/Rscripts/120704-sortDataFrame.R')
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')

#Draw a heatmap. First normalise the raw counts by library size
dat = read.delim('140213_normalisedLog_CPM_facs.txt')
datM = as.matrix(dat)
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

##################### Extract the top 2000 most variable probes ###########################
colour = c("forestgreen", "forestgreen", "forestgreen", "blue", "blue", "blue")
colour = as.factor(colour)

# Extract median absolute deviation and take the top 500
madData = apply(datM, 1, mad)

madDataSort = as.data.frame(madData)
madDataSort$probe = row.names(madDataSort)
colnames(madDataSort) = c('value', 'probe')
madDataSort = sort.dataframe(madDataSort, 1, highFirst=T)
top500 = row.names(madDataSort[c(1:1500),])
topNorm = datM[top500,]
# name = paste(dmCut$sample, dat$Replicate)

###################### Make heat map with Veerhaak subtype #####################
heatmap.2(topNorm, cexRow=1.5, main="", #scale='row',
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", 
          labRow=NA, xlab="Samples", labCol=colnames(topNorm),
          offsetRow=c(1,1), margins=c(8,4))
