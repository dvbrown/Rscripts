library(DEXSeq)
library(multicore)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

setwd('~/Documents/RNAdata/danBatch1/dexSeq_count/')
samples = c('long1', 'long2', 'long3', 'short1', 'short2', 'short3')
dm = read.csv('../bowtieGem/revHTSeq/designMatrix.csv')
dm = dm[,c(1,2,4,5)]
condition = as.factor(c('long', 'long', 'long', 'short', 'short', 'short'))
#design matrix
data = read.HTSeqCounts(c('GIC_011.countExons.txt', 'GIC_020.countExons.txt', 'GIC_034.countExons.txt', 'GIC_035.countExons.txt',
                    'GIC_039.countExons.txt', 'GIC_041.countExons.txt'), condition,
                    flattenedfile='~/Documents/RNAdata/pilotRNAseq/121121_TopHatNovelAlignment/121203_countExons/exonAnnotation.gtf')


head(fData(data))
data = estimateSizeFactors(data)
sizeFactors(data) #scale for library size
data = estimateDispersions(data, minCount=20) #there is no replicates so will have to manually set one.
# fData(data)$dispersion = 0.1
# fData(data)$dispFitted = 0.1
data = fitDispersionFunction(data)
head(fData(data)$dispBeforeSharing)
fData(data)$testable = ifelse((rowSums(counts(data) > 50)), TRUE, FALSE) #only test exon usage for genes with more than 50 counts
data@dispFitCoefs
head(fData(data)$dispFitted)
plotDispEsts(data)

condition = as.factor(c('long', 'long', 'long', 'short', 'short', 'short'))
data@phenoData@data$condition

#data = testForDEU(data)#, formula0=count ~ sample + exon + data[[2]], formula1=count ~ sample + exon + data[[2]] *I (exon==exonID))
data = testForDEUTRT(data)
dataNorm = estimatelog2FoldChanges(data)
result = DEUresultTable(data)
table(result$padjust != NA)

plot(result$meanBase, result[, "log2fold(negative/positive)"], log = "x",
     col = ifelse(res1$padjust < 0.1, "red", "black"), ylim = c(-4,4), main = "CD133 MvsA")

keep = rowSums(counts(dataNorm) > 10)