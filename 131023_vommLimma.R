library(edgeR)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/')
source('~/Documents/Rscripts/annotateEnsembIDs.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')
files = list.files(pattern='*.txt')

dm = read.csv('designMatrix.csv')
f = lapply(files, read.delim, header=FALSE)
df1 = cbind(f[[1]],f[[2]],f[[3]],f[[4]],f[[5]],f[[6]])
df = df1[,c(2,4,6,8,10,12)]
row.names(df) = df1[,1]
colnames(df) = c('GIC_011', 'GIC_020', 'GIC_034', 'GIC_035', 'GIC_039', 'GIC_041')
labels = c('#011', '#020', '#034', '#035', '#039', '#041')

noFeatures = tail(df)
noCount = rownames(df) %in% c("no_feature","ambiguous","too_low_aQual",
                              "not_aligned","alignment_not_unique")
df = df[!noCount,]
totalCount = colSums(df)

#Build EdgeR objects
condition = c('long', 'long', 'long', 'short', 'short', 'short')
counts = DGEList(counts=df, group=condition)

# I remove genes with less than 0.5 cpm in 3 samples. For a library size of 20M this is 10 reads.
cpms = cpm(counts)
keep = rowSums(cpms >0.5) >=3
counts = counts[keep,]

############################## The Voom Limma Part ##############################################
#apply normalisation
counts = calcNormFactors(counts)
#Use voom to convert the read counts to log2-cpm, with associated weights, ready for linear modelling:
design = model.matrix(~ libPrep + age + group, dm)
v <- voom(counts, design, plot=TRUE)

#Make a MDS plot to view differences
plotMDS(v,top=500, col = c("darkgreen","blue")[factor(condition)], labels=row.names(counts$samples), 
        gene.selection="pairwise", main='MDS plot RNA-seq batch1', cex=1.25)

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
