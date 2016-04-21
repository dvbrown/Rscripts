library(GSVA)
library(gplots)

##################### IO ###########################
setwd('~/Documents/public-datasets/GPSC_subgroups/bhat2013/analysis/')

data = read.delim('140826_probeAveraged.txt', row.names=1)
dataM = as.matrix(data)
dm = read.delim('../designMatrix.txt')[c(1:17),]
# Just the untreated controls

cd133Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd133Cutoff.txt", row.names=1)
cd44Sig = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/140527_cd44Cutoff.txt", row.names=1)
cd15 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15/140528_cd15Cutoff.txt", row.names=1)
aldh1 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/140528_ALDH1Cutoff.txt", row.names=1)
itag6 = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6//140528_ITGA6Cutoff.txt", row.names=1)
l1cam = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/140528_L1CAMCutoff.txt", row.names=1)

bigSigs = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), "CD15" = row.names(cd15),
               "ALDH1"=row.names(aldh1), "ITGA6"=row.names(itag6), "L1CAM"=row.names(l1cam))
rm(cd133Sig, cd44Sig, cd15, aldh1, itag6, l1cam)

##################### GSVA ###########################
bigResult = gsva(dataM, bigSigs,  rnaseq=F, verbose=T, parallel.sz=1)
bigResult = t(bigResult$es.obs)

dm$colour = "black"
dm$colour[dm$Subtype %in% 'Proneural'] = 'purple'
dm$colour[dm$Subtype %in% 'Mesenchymal'] = 'red'
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

##################### heatMaps ###########################

heatmap.2(t(bigResult), cexRow=1.2, main="Enrichment of FACS marker signatures \n in Molecular Subtype", scale="none",
          Rowv=NULL, Colv=TRUE, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", 
          ColSideColors=as.character(dm$colour), xlab="Bhat et al 2013 samples", labCol=row.names(bigResult), labRow=colnames(bigResult), 
          offsetRow=c(1,1), margins=c(8,7))

smallResult = bigResult[,c(1:3)]
heatmap.2(t(smallResult), cexRow=1.2, main="Enrichment of FACS marker signatures \n in Molecular Subtype", scale="none",
          Rowv=NULL, Colv=TRUE, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", 
          ColSideColors=as.character(dm$colour), xlab="Bhat et al 2013 samples", labCol=row.names(smallResult), labRow=colnames(smallResult), 
          offsetRow=c(1,1), margins=c(8,7))

# Do single genes
markers = as.matrix(data[c('CD44', 'FUT4', 'PROM1'),])
heatmap.2(markers, cexRow=1.2, main="Enrichment of FACS marker mRNAs \n in Molecular Subtype", scale="none",
          Rowv=NULL, Colv=TRUE, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", 
          ColSideColors=as.character(dm$colour), xlab="Bhat et al 2013 samples", labCol=colnames(markers), labRow=c('CD44', 'CD15', 'CD133'), 
          offsetRow=c(1,1), margins=c(8,7))