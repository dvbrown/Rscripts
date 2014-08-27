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

setwd('~/Documents/public-datasets/GPSC_subgroups/bhat2013/GSE49009_RAW/')
list.files()

rawData = ReadAffy()
dm = readTargets('../designMatrix.txt')
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

par(mfrow=(c(2,1)))
# boxData = exprs(rawData)
# colnames(boxData) = dm$cellLine
boxplot(rawData, col=rainbow(21), main='Bhat unNormalised data', las=1, mar=c(5,5,5,2))
data = rma(rawData, verbose=T)
norm = exprs(data)
colnames(norm) = dm$cellLine
boxplot(norm, col=rainbow(21), main='Bhat normalised data', las=2, mar=c(5,5,5,2))
par(mfrow=(c(1,1)))

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

g2 = ggplot(data=pcaDf, aes(x=V1, y=V2, color=Source)) + 
    geom_point(shape=19) + #geom_smooth(method=lm, colour='red') +
    xlab("PC1") + ylab("PC2") + # Set axis labels
    ggtitle("Principal component analysis Bhat 2013") +  # Set title
    theme_bw(base_size=18)

g2 + geom_text(data = pcaDf, aes(x=V1, y=V2,
                                label=dd_text[1:17], size=0.2), colour="black")

##################### Make a heatmap ###########################
dm$colour = "black"
dm$colour[dm$Subtype %in% 'Proneural'] = 'purple'
dm$colour[dm$Subtype %in% 'Mesenchymal'] = 'red'

# Extract median absolute deviation and take the top 500
madData = apply(norm, 1, mad)

madDataSort = as.data.frame(madData)
madDataSort$probe = row.names(madDataSort)
colnames(madDataSort) = c('value', 'probe')
madDataSort = sort.dataframe(madDataSort, 1, highFirst=T)
top500 = row.names(madDataSort[c(1:500),])
topNorm = norm[top500,]
topNorm = sort.dataframe(topNorm, 1, highFirst=T)

heatmap.2(topNorm[,c(1:17)], cexRow=0.8, main="Enrichment of FACS marker signatures \n in Molecular Subtype", scale="row",
          Rowv=NULL, Colv=TRUE, keysize=1, trace="none", col=myPalette, density.info="none", dendrogram="column", 
          ColSideColors=as.character(dm$colour)[1:17], labRow=NA, xlab="Aglent samples", labCol=colnames(topNorm), 
          offsetRow=c(1,1), margins=c(10,4))

##################### Export data ###########################
write.table(norm, '../analysis/140826_rmaNormalised.txt', sep='\t')

##################### Differential expression testing ###########################
dataUntreated = norm[,c(1:17)]

f = paste(dm$Subtype)[c(1:17)]
f = factor(f)
design = model.matrix(~0+f)
colnames(design) = levels(f)

fit = lmFit(dataUntreated, design)
fit$design

cont.matrix = makeContrasts(mesVSpn="Mesenchymal-Proneural", levels=design)
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
summariseData = summariseData[row.names(summariseData) != 'NA' ,c(1:17)]
write.table(summariseData, '../analysis/140826_probeAveraged.txt', sep='\t')
