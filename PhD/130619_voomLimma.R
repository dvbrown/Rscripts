library(limma)
library(edgeR)
library(biomaRt)
setwd('~/Documents/RNAdata/RNAseqAnalysis/')

#retrieve EMSEMBL gene ID mappings from biomaRt
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
#add the annotations to the dataset
data = read.delim('121030_readCountsGenomeRanges.txt', row.names=1)
libSize = colSums(data)
groups=c('negative','negative','negative','positive', 'positive', 
         'positive','unsorted', 'unsorted', 'unsorted','test','test')
annotation = annotateIds(data)

#build the DGE object
diffGeneExpress <- DGEList(counts=data, lib.size=libSize, group=groups)
# filter for genes with at least 1 count per million  in at least 3 samples
isexpr <- rowSums(cpm(diffGeneExpress)>1) >= 3
filData = diffGeneExpress[isexpr,]

#apply normalisation
filData = calcNormFactors(filData)
#Use voom to convert the read counts to log2-cpm, with associated weights, ready for linear modelling:
design <- model.matrix(~filData$samples$group)
v <- voom(filData,design,plot=TRUE)

#Make a MDS plot to view differences
plotMDS(v,top=50,labels=filData$samples$group, gene.selection="common", main='MDS RNA-seq \nclone 035 pilot')

#differential exppression test as for limma
fit <- lmFit(v,design)
fit <- eBayes(fit)
#inputting the gene list doesn't generate the right annotation ENSG mappings
posVSneg = topTable(fit,coef=2,number=Inf,sort.by="p", resort.by='logFC')
posVSneg = merge(posVSneg, annotation, by.x='ID',by.y='ensembl_gene_id')
head(summary(decideTests(fit)))
write.table(posVSneg, '~/Documents/RNAdata/RNAseqAnalysis/121105_trimmomaticReads/mergedBam/121107_mergeSortTopHatAlignIndex/130619_voomLimma/130620_voomLimma.txt',sep='\t',row.names=F)

#filter for known genes, significant genes and differentially expressed genes. Use a generic subset idea, take all columns
filPosVSneg = posVSneg[(abs(posVSneg$logFC) > 1) & (posVSneg$adj.P.Val < 0.05),]
write.table(filPosVSneg, '~/Documents/RNAdata/RNAseqAnalysis/121105_trimmomaticReads/mergedBam/121107_mergeSortTopHatAlignIndex/130619_voomLimma/130620_voomLimmaFilter.txt',sep='\t',row.names=F)
     
plotMA(fit,array=2,main="MA plot clone 035 pilot voom limma", xlab='logCounts', ylab='logFC')