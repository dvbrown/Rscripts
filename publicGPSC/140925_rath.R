library(limma)
library(affy)
library(gplots)
library(ggplot2)
source('~/Documents/Rscripts/120704-sortDataFrame.R')
source('~/Documents/Rscripts/multiplot.R')

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
setwd('~/Documents/public-datasets/GPSC_subgroups/rath2012_coCultureAstro/GSE37120_RAW/')
list.files()

rawData = ReadAffy()
dm = readTargets('../design.txt')
dm$cellLine = paste(dm[,2], dm$suorce, dm$subpopulation)
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

par(mfrow=(c(2,1)))
boxData = exprs(rawData)
colnames(boxData) = dm$cellLine
boxplot(rawData, col=rainbow(21), main='Rath unNormalised data', las=1, mar=c(5,5,5,2))
data = rma(rawData, verbose=T)
norm = exprs(data)
colnames(norm) = dm$cellLine
boxplot(norm, col=rainbow(21), main='Rath normalised data', las=2, mar=c(5,5,5,2))
par(mfrow=(c(1,1)))

##################### Remove non-informatives entries ###########################
dm$sample = paste(dm$Source, dm$Subpopulation)
dmCut = dm[c(5,6,11,12,17,18,23,24),]
norm = norm[,c(5,6,11,12,17,18,23,24)]

##################### Render the prinical components ###########################
pca = princomp(norm)
pcaDf = as.data.frame(cbind(pca$loadings[,1], pca$loadings[,2]))
pcaDf = cbind(pcaDf, dmCut)
dd_text = dm$sample

g = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Source)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Ruth 2012") +  # Set title
    theme_bw(base_size=18)
g
g2 = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Subpopulation)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Ruth 2012") +  # Set title
    theme_bw(base_size=18)
g2
g3 = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Replicate)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Ruth 2012") +  # Set title
    theme_bw(base_size=18)
g3
multiplot(g, g2,g3, cols=2)

##################### Make a heatmap ###########################
dmCut$colour = c("blue", "blue", "red", "red", "blue", "blue", "red", "red")
dmCut$colour = as.factor(dmCut$colour)

# Extract median absolute deviation and take the top 500
madData = apply(norm, 1, mad)

madDataSort = as.data.frame(madData)
madDataSort$probe = row.names(madDataSort)
colnames(madDataSort) = c('value', 'probe')
madDataSort = sort.dataframe(madDataSort, 1, highFirst=T)
top500 = row.names(madDataSort[c(1:500),])
topNorm = norm[top500,]
# topNorm = sort.dataframe(topNorm, 1, highFirst=T)
name = paste(dmCut$sample, dmCut$Replicate)

heatmap.2(topNorm, cexRow=0.8, main="Gene expression profiles CD133 sorted GPSCs", scale="row",
          Colv=as.factor(dmCut$colour), keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="row", 
          ColSideColors=as.character(dmCut$colour), labRow=NA,labCol=name, 
          offsetRow=c(1,1), margins=c(10,4))

heatmap.2(topNorm, cexRow=0.8, main="Gene expression profiles CD133 sorted GPSCs", scale="row",
          keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="both", 
          ColSideColors=as.character(dmCut$colour), labRow=NA, labCol=name, 
          offsetRow=c(1,1), margins=c(15,4))

# Take the probe averages
summariseData = avereps(data, ID=genelist$GeneSymbol)
# Remove NAs
summariseData = summariseData[row.names(summariseData) != 'NA' ,]

##################### Export data ###########################
write.table(norm, '../analysis/140911_rmaNormalised.txt', sep='\t')
write.table(summariseData, '../analysis/140911_ruthGEM.txt', sep='\t')
write.table(topNorm, '../analysis/140911_top500.txt', sep='\t')