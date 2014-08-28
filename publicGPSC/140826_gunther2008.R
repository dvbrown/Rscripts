library(limma)
library(affy)
library(gplots)
library(ggplot2)
source('~/Documents/Rscripts/120704-sortDataFrame.R')
##################### Build the genelist ###########################

library(annotate)
library(hgu133a2.db)
geneIDs = ls(hgu133a2ENTREZID)

#get gene ID numbers from the annptation package allowing for multiple probes to match mulitple genes
geneSymbols <- as.character(unlist(lapply(mget(geneIDs,env=hgu133a2SYMBOL),
                                          function (symbol) { return(paste(symbol,collapse="; ")) } )))
geneNames <- as.character(unlist(lapply(mget(geneIDs,env=hgu133a2GENENAME),
                                        function (name) { return(paste(name,collapse="; ")) } )))
unigene <- as.character(unlist(lapply(mget(geneIDs,env=hgu133a2UNIGENE),
                                      function (unigeneID) { return(paste(unigeneID,collapse="; ")) } )))

#strip the Hs from the start of unigene reference
unigene <- gsub("Hs\\.","",unigene)

#read the gene annotations into a dataframe for use in the topTable function of limma
genelist <- data.frame(GeneID=geneIDs,GeneSymbol=geneSymbols,GeneName=geneNames)
rm(geneSymbols, geneNames, unigene)

##################### IO ###########################

setwd('~/Documents/public-datasets/GPSC_subgroups/Gunther2008-Oncogene/GSE8049_RAW/')
list.files()

rawData = ReadAffy()

dm = readTargets('../designMatrix.txt')
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

par(mfrow=(c(2,1)))
# boxData = exprs(rawData)
colnames(boxData) = dm$cellLine
boxplot(rawData, col=rainbow(21), main='Gunther unNormalised data', las=1, mar=c(5,5,5,2))
# data = rma(rawData, verbose=T)
# norm = exprs(data)
colnames(norm) = dm$cellLine
boxplot(norm, col=rainbow(21), main='Gunther normalised data', las=2, mar=c(5,5,5,2))
par(mfrow=(c(1,1)))

rm(boxData)
##################### Render the prinical components ###########################
pca = princomp(norm)
pcaDf = as.data.frame(cbind(pca$loadings[,1], pca$loadings[,2]))
pcaDf = cbind(pcaDf, dm)[c(1:17),]
dd_text = dm$cellLine

g = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Subtype)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Bhat 2013") +  # Set title
    theme_bw(base_size=18)

g + geom_text(data = pcaDf, aes(x=V1, y=V2,
                                label=dd_text[1:17], size=0.2), colour="black")
