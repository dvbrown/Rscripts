library(limma)
library(affy)
library(gplots)
library(ggplot2)
source('~/Documents/Rscripts/120704-sortDataFrame.R')
##################### Build the genelist ###########################

library(annotate)
library(hgu133plus2.db)
geneIDs = ls(hgu133plus2ENTREZID)

#get gene ID numbers from the annptation package allowing for multiple probes to match mulitple genes
geneSymbols <- as.character(unlist(lapply(mget(geneIDs,env=hgu133plus2SYMBOL),
                                          function (symbol) { return(paste(symbol,collapse="; ")) } )))
geneNames <- as.character(unlist(lapply(mget(geneIDs,env=hgu133plus2GENENAME),
                                        function (name) { return(paste(name,collapse="; ")) } )))
unigene <- as.character(unlist(lapply(mget(geneIDs,env=hgu133plus2UNIGENE),
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

g = ggplot(data=pcaDf, aes(x=V1, y=V2, color=subtype)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Gunther 2008") +  # Set title
    theme_bw(base_size=18)

g + geom_text(data = pcaDf, aes(x=V1, y=V2,
                                label=dd_text[1:17], size=0.2), colour="black")

##################### Make a heatmap ###########################
dm$colour = "black"
dm$colour[dm$subtype %in% 'cluster1'] = 'purple'
dm$colour[dm$subtype %in% 'cluster2'] = 'red'

# Extract median absolute deviation and take the top 500
madData = apply(norm, 1, mad)

madDataSort = as.data.frame(madData)
madDataSort$probe = row.names(madDataSort)
colnames(madDataSort) = c('value', 'probe')
madDataSort = sort.dataframe(madDataSort, 1, highFirst=T)
top500 = row.names(madDataSort[c(1:500),])
topNorm = norm[top500,]
topNorm = sort.dataframe(topNorm, 1, highFirst=T)

heatmap.2(topNorm, cexRow=0.8, main="Clustering of Gunther et al microarray", scale="row",
          Rowv=NULL, Colv=TRUE, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", 
          ColSideColors=as.character(dm$colour), labRow=NA, xlab="Samples", labCol=colnames(topNorm), 
          offsetRow=c(1,1), margins=c(7,4))

##################### Export data ###########################
write.table(norm, '../analysis/140826_rmaNormalised.txt', sep='\t')

##################### Differential expression testing ###########################

f = paste(dm$subtype)
f = factor(f)
design = model.matrix(~0+f)
colnames(design) = levels(f)

fit = lmFit(norm, design)
fit$design

cont.matrix = makeContrasts(mesVSpn="cluster2-cluster1", levels=design)
fit2  = contrasts.fit(fit, cont.matrix)
fit2  = eBayes(fit2)

#write the output to a table of differentially expressed genes. Change this value to suit
result = topTable(fit2, coef=1, number=22277, genelist=genelist,  
                  adjust='BH', sort.by='logFC', lfc=0)
result = sort.dataframe(result, 4, highFirst=T)

head(result,100)
tail(result,100)

results <- decideTests(fit2)
summary(results)

write.table(result, '../analysis/140826_mesVSpnDE.txt', sep='\t')

# Tkae the probe averages
summariseData = avereps(data, ID=genelist$GeneSymbol)
colnames(summariseData) = dm$cellLine

# Remove NAs
summariseData = summariseData[row.names(summariseData) != 'NA' ,]
write.table(summariseData, '../analysis/140826_probeAveraged.txt', sep='\t')
