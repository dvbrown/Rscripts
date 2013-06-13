#Enrichment score. For example, 10% of user’s genes are kinases versus 1% of genes in human genome 
#(this is population background) are kinases. Thus, the fold enrichment is tenfold. 
#Fold enrichment along with EASE score could rank the enriched terms in a more comprehensive way
####################################################################################################
setwd('~/Documents/CREB/ChIPseqENCODE/DAVID/')
data = read.delim('GeneOntogeny.txt')
data = data[,c(1,2,3,4,7,8,9,10,11)]
colnames(data) = c('category', "term", 'count', 'percentTotal', 'total', 'popHits', 'popTotal', 'enrichment', 'FDRcorrected')
data1 = gsub('GO:.......~', '', data$term)
data$term = data1
data$category = as.character(data$category)
#sort the data remeber to include the comma for columns!
sortData = data[order(data$category),]

#the filters
filterData = subset(sortData, sortData$FDRcorrected <= 0.1)
filterData = subset(filterData, filterData$enrichment >= 2)

goTerm = filterData[1:43,]
goTerm = goTerm[order(goTerm$enrichment),]
goTerm[42,2] = 'RNA polymerase II transcription factor activity'
goTerm = goTerm[c(-1,-15,-16,-18,-19,-20,-30,-35,-39),]
spIR = filterData[53:68,]
spIR = spIR[order(spIR$enrichment),]

par(mar=c(4,8.5,5,2))
barplot(spIR$enrichment, names.arg=spIR$term, horiz=T, las=2, cex.axis=0.8, main='Enriched processes CREB ENCODE ChIP',
        xlab='Fold enrichment', col=rainbow(16))

par(mar=c(4,13,5,2))
barplot(goTerm$enrichment, names.arg=goTerm$term, horiz=T,las=2, cex.axis=0.8, main='Enriched molecular functions CREB ENCODE ChIP',
        xlab='Fold enrichment', col=rainbow(12), cex=0.7)

par(mfrow=c(2,1))
#Plot the transcription factor data
tfs = read.delim('transcription factorBinding.txt')
tfs = tfs[,c(2,3,4,8,10,13)]
tfs = subset(tfs, tfs$FDR <= 0.05)
tfs = subset(tfs, tfs$Fold.Enrichment >= 1.2)
tfs = tfs[order(tfs$Fold.Enrichment, decreasing=F),]

#Plot tissue expression profiles
tissue = read.delim('tissueExpression.txt')
tissue = tissue[,c(2,3,4,8,10,13)]
tissue = subset(tissue, tissue$FDR <= 0.05)
tissue = tissue[order(tissue$Fold.Enrichment, decreasing=F),]

barplot(tfs$Fold.Enrichment, names.arg=tfs$Term, horiz=T,las=2, cex.axis=0.8, main='Enriched transcription \nfactors CREB ChIP',
        xlab='Fold enrichment', col=rainbow(12), cex=0.66)
barplot(tissue$Fold.Enrichment, names.arg=tissue$Term, horiz=T,las=2, cex.axis=0.8, main='CREB ChIP tissue \nexpression profile',
        xlab='Fold enrichment', col=rainbow(12), cex=0.6)