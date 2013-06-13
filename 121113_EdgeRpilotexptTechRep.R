library(edgeR)
setwd('~/Documents/RNAdata/RNAseqAnalysis/')
data = read.delim('121030_readCountsGenomeRanges.txt', row.names=1)
data = data[,1:6]
#sum up the technical replicates for a single table of counts
CD133neg = rowSums(data[,c(1:3)])
CD133pos = rowSums(data[,c(4:6)])
CD133unsort = rowSums(data[,c(7:9)])
CD133test = rowSums(data[,c(10,11)])
sumData = cbind(CD133neg, CD133pos, CD133unsort, CD133test)
dataDesign =c('negative','positive', 'unsorted','test')
dge = DGEList(counts=sumData, group=dataDesign)
normFact = calcNormFactors(dge)
#MA plots of pre and post normalised data
maPlot(normFact$counts[,1], normFact$counts[,2], normalize=TRUE, pch=19, cex=0.4, ylim=c(-8, 8))
grid(col = "blue")
abline(h = log2(normFact$samples$norm.factors[2]/normFact$samples$norm.factors[1]), col="red", lwd=4)
eff.libsize <- normFact$samples$lib.size * normFact$samples$norm.factors
maPlot(normFact$counts[, 1]/eff.libsize[1], normFact$counts[, 2]/eff.libsize[2], normalize = FALSE, pch = 19, cex = 0.4, ylim = c(-8, 8))
grid(col = "blue")

#classical one factor analysis. As technical replicates can't fit dispersion. Use the recommended value from documentation
et = exactTest(normFact, pair=c('negative','positive'), dispersion=0.1)
deGenes = topTags(et, 20000, sort.by='p.value')
finalData = deGenes$table

#Try the general linear model
replicates = c('neg1', 'neg2', 'neg3', 'pos1', 'pos2', 'pos3')
samples = c('CD133low', 'CD133low', 'CD133low', 'CD133high', 'CD133high', 'CD133high')
design = model.matrix(~replicates+samples) #design matrix
dataDesign =c('negative','negative','negative','positive', 'positive','positive')
dge = DGEList(counts=data, group=dataDesign)
dge = calcNormFactors(dge)
disp = estimateGLMCommonDisp(dge, design)
disp = estimateGLMTrendedDisp(dge, design)
fit = glmFit(dge, design, dispersion=0.1 )
#Design matrix not of full rank. There is not a second experimental factor to examine or even a batch number.

#Have to adjust this function for edgeR
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
annotateIds = function(geneList)  { #Given a gene list return the ENSEMBL information and ID mappings
  ensembl_genes = row.names(geneList)
  list=getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "external_gene_id", "entrezgene", "description"),
    values= ensembl_genes,
    mart= mart)
  return (list)
}
geneList = as.data.frame(deGenes)
annResult = annotateIds(deGenes)
finalResult = merge.data.frame(deGenes, annResult, by.x='row.names', by.y='ensembl_gene_id')
write.table(finalResult, './121105_trimmomaticReads/mergedBam/121107_mergeSortTopHatAlignIndex/121121_edgeRanalysisNoReps/121122_edgeRclassical.txt', sep='\t', row.names=T)