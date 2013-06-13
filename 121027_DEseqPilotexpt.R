library(DESeq)
library(ggplot2)
library(Heatplus)
setwd('~/Documents/RNAdata/RNAseqAnalysis/')
data = read.delim('121030_readCountsGenomeRanges.txt', row.names=1)
plot(data$CD133n_rep2, data$CD133n_rep3, xlim=c(0,25000), ylim=c(0,25000), cex=0.25, type='p', ylab='CD133 low rep1', xlab='CD133 low rep2', main='Read Count correlation replicates')
lines(lowess(data$CD133n_rep1, data$CD133n_rep3), col="blue") # lowess line (x,y)
plot(data$CD133p_rep1, data$CD133p_rep3, xlim=c(0,25000), ylim=c(0,25000), cex=0.25, type='p', ylab='CD133 high rep2', xlab='CD133 high rep3', main='Read Count correlation replicates')
lines(lowess(data$CD133p_rep2, data$CD133p_rep3), col="blue") # lowess line (x,y)
plot(data$CD133x_rep1, data$CD133x_rep2, xlim=c(0,25000), ylim=c(0,25000), cex=0.25, type='p', ylab='CD133 unsorted rep2', xlab='CD133 unsorted rep3', main='Read Count correlation replicates')
lines(lowess(data$CD133x_rep2, data$CD133x_rep1), col="blue") # lowess line (x,y)

#Design matrix
dataDesign = data.frame(row.names=colnames(data)
                        ,condition=c('negative','negative','negative','positive','positive', 'positive',
                                     'unsorted', 'unsorted', 'unsorted', 'test','test'),
                        libType=c('paired-end','paired-end','paired-end','paired-end','paired-end','paired-end',
                                  'paired-end','paired-end','paired-end','paired-end','paired-end'))

conds = dataDesign$condition
#create a count data set object
countDataSet = newCountDataSet(data, conds)

#Get the library size scaled data
cds = estimateSizeFactors(countDataSet)
normaliseData = counts(cds, normalized=T)
#fiddle with the dispersion calculation if you don't get anything significant
cds = estimateDispersions(cds, method='pooled')

#A function to draw a dispersion vs normalised read count
plotDispEsts <- function( countDataSet ) {
  plot(
    rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds)$perGeneDispEsts,
    pch = '.', log="xy", main='Fitted disperson parameter', ylab='Dispersion estimate per gene', xlab='normalised read count' )
  xg <- 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}

plotDE <- function( res ) { #Draw the MA plot
  plot(
    res$baseMean,
    res$log2FoldChange,
    log="x", pch=20, cex=.3,
    col = ifelse( res$padj < .05, "red", "black" ), main='MA plot CD133 low vs CD133 high', ylab='log2FC', xlab='read count' )
}

plotDispEsts(cds) #Plot the dispersion

result = nbinomTest(cds, 'negative', 'positive') #the differential expression test
result.unsort.positive = nbinomTest(cds, 'unsorted', 'positive') #the differential expression test
result.unsort.negative = nbinomTest(cds, 'unsorted', 'negative') #the differential expression test

plotDE(result) #MA plot
hist(result$pval, breaks=100, col="skyblue", border="slateblue", main="Significance of differential expression tests", xlab='FDR adjusted p-value') #histogram the results of the differential expression test

sigResult = result[result$padj < 0.1,] #subset the most significant genes
head( sigResult[ order(sigResult$pval), ] ) #list most significant genes

ncuRep <- counts( cds, normalized=TRUE )[ , conditions(cds)=="negative" ] #extract replicate normalised counts
plot( rowMeans(ncuRep), log2( ncuRep[,2] / ncuRep[,1] ), pch=".", log="x" ) #plot the normalised replicate counts

require(ggplot2)
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
finalResult$threshold = as.factor(abs(finalResult$log2FoldChange) > 2 & finalResult$padj < 0.05)

##Construct the volcano plot object. First get the final result by running annotate ENSEMBL IDs script with genelist
g = ggplot(data=finalResult, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  opts(legend.position = "none" , title="Differential expression CD133 low vs CD133 high glioma initiating cells" , plot.title=theme_text(size=16)) +
  #xlim(c(-5, 5)) + ylim(c(0, 50)) +
  xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value")
g

#subset gene names for only significant genes
dd_text = finalResult[(abs(finalResult$log2FoldChange) > 4.33) & (finalResult$padj < 0.05),]
#add text to volcano
g + geom_text(data = dd_text, aes(x=log2FoldChange, y=-log10(padj),
                  label=external_gene_id, size=0.25), colour="black")

#Draw a heatmap. First normalise the raw counts by library size
data = read.delim('121029_sharingModeFitOnly/CD133lowVSCD133hi_cutoffs/121029_DE_CD133lowCD133hi_ApValLFC_Cutoff1-5.txt')
row.names(data) = data$id
rawData = read.delim('121030_readCountsGenomeRanges.txt')
sigGenes = as.character(data$id) #use to index the raw data

dataMatrix = data.frame(row.names=colnames(rawData)
                        ,condition=c('CD133low_Rep1','CD133low_Rep2','CD133low_Rep3','CD133high_Rep1','CD133high_Rep1', 'CD133high_Rep1',
                                     'CD133unsort_rep1', 'CD133unsort_rep2', 'CD133unsort_rep3', 'test1','test2'),
                        libType=c('paired-end','paired-end','paired-end','paired-end','paired-end','paired-end',
                                  'paired-end','paired-end','paired-end','paired-end','paired-end'))
conds = dataMatrix$condition
countDataSet = newCountDataSet(rawData, conds)
#Get the library size scaled data
cds = estimateSizeFactors(countDataSet)
normaliseData = counts(cds, normalized=T) #use normalised data for the heat map
sigData = normaliseData[sigGenes,]
geneNames = data$external_gene_id
row.names(sigData) = geneNames
#Now call the heat map
corrdist = function(x) as.dist(1-cor(t(x)))
heatmap(sigData, margins=c(7,5))
