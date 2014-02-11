library(DEXSeq)
library(multicore)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

setwd('~/Documents/RNAdata/danBatch1/dexSeq_count/')
samples = c('long1', 'long2', 'long3', 'short1', 'short2', 'short3')
# Read in design matrix
dm = read.csv('../bowtieGem/revHTSeq/designMatrix.csv')
dm = dm[,c(1,2,4,5)]
condition = as.factor(c('long', 'long', 'long', 'short', 'short', 'short'))

# Read the output of DEXSeq count into the appropriate object
data = read.HTSeqCounts(c('GIC_011.countExons.txt', 'GIC_020.countExons.txt', 'GIC_034.countExons.txt', 'GIC_035.countExons.txt',
                    'GIC_039.countExons.txt', 'GIC_041.countExons.txt'), condition,
                    flattenedfile='~/Documents/RNAdata/pilotRNAseq/121121_TopHatNovelAlignment/121203_countExons/exonAnnotation.gtf')


head(fData(data))
data = estimateSizeFactors(data)
# Scale for library size
sizeFactors(data) 

# Estimate dispersion
data = estimateDispersions(data, minCount=20)
data = fitDispersionFunction(data)
head(fData(data)$dispBeforeSharing)
fData(data)$testable = ifelse((rowSums(counts(data) > 50)), TRUE, FALSE) #only test exon usage for genes with more than 50 counts
data@dispFitCoefs
head(fData(data)$dispFitted)
# Can't find the fuction -> plotDispEsts(data)

data@phenoData@data$condition

# Test for differential expression using the TRT method
data = testForDEUTRT(data)
data= estimatelog2FoldChanges(data)
result = DEUresultTable(data)
table(result$padjust != NA)

# Sort the dataframe for p-value
resultSorted = sort.dataframe(result, 5, highFirst=FALSE)
colnames(resultSorted) = c("geneID","exonID","dispersion","pvalue","padjust",'meanBase', 'lfc_short'  )

# Subset the dataframe for interesting genes
sigGenes = resultSorted[(resultSorted$padjust < 0.1) & (abs(resultSorted$lfc_short) >= 1), ]

# This plot function doesn't work well
plot(resultSorted$meanBase, resultSorted$lfc_short, log = "x",
     col = ifelse(resultSorted$padjust < 0.1, "red", "black"), ylim = c(-4,4), main = "CD133 MvsA")

# Write the results out to file

write.table(resultSorted, '131114_dexSeqResultSorted.txt', sep ='\t')
write.table(sigGenes, '131114_dexSeqSigGenes.txt', sep ='\t')

plotDEXSeq(data, 'ENSG00000206503+ENSG00000227766', legend=TRUE, color=c('darkgreen', 'blue'), 
           names=TRUE, expression=TRUE, main='HLA-A + HLA complex group 4')

plotDEXSeq(data, 'ENSG00000214940+ENSG00000205746+ENSG00000233024+ENSG00000254681', legend=TRUE, color=c('darkgreen', 'blue'), 
           names=TRUE, expression=TRUE, main='Various psuedogenes')

plotDEXSeq(data, 'ENSG00000221988+ENSG00000241404+ENSG00000204314+ENSG00000204310+ENSG00000258388', legend=TRUE, color=c('darkgreen', 'blue'), 
           names=TRUE, expression=TRUE, main='PPT2 + EGFL8 + PRRT1 + AGPAT1 + PPT2-EGFL8 readthrough')

plotDEXSeq(data, 'ENSG00000126822+ENSG00000070182', legend=TRUE, color=c('darkgreen', 'blue'), 
           names=TRUE, expression=TRUE, norCounts=F, main='PLEKHG3 + SPTB')

plotDEXSeq(data, 'ENSG00000227199+ENSG00000226367+ENSG00000214188+ENSG00000004866', legend=TRUE, color=c('darkgreen', 'blue'), 
           names=TRUE, expression=TRUE, norCounts=F, main='ST7 antisense RNA 1 + 2 + ST7 overlapping transcript 4 + ST7')

plotDEXSeq(data, 'ENSG00000129675', legend=TRUE, color=c('darkgreen', 'blue'), 
           names=TRUE, expression=TRUE, norCounts=F, main='Rac/Cdc42 guanine nucleotide exchange factor (GEF) 6')

plotDEXSeq(data, 'ENSG00000163617+ENSG00000151576+ENSG00000184307', legend=TRUE, color=c('darkgreen', 'blue'), 
           names=TRUE, expression=TRUE, norCounts=F, main='KIAA1407 + QTRTD1 + ZDHHC23')
