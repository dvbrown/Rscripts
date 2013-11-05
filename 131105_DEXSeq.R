library(DEXSeq)
library(multicore)

setwd('~/Documents/RNAdata/RNAseqAnalysis/121121_TopHatNovelAlignment/121203_countExons/')
samples = as.character(c('negative', 'positive'))
#design matrix
data = read.HTSeqCounts(c('s_4_CD133n_A.bam.sorted.mergeBam.addReadGroupHeader.reorderSam.markDuplicates.sortReadName.countExons.txt', 
                        's_4_CD133p_A.bam.sorted.mergeBam.addReadGroupHeader.reorderSam.markDuplicates.sortReadName.countExons.txt'), 
                    samples, flattenedfile='exonAnnotation.gtf')
head(fData(data))
data = estimateSizeFactors(data)
sizeFactors(data) #scale for library size
data = estimateDispersions(data) #there is no replicates so will have to manually set one.
fData(data)$dispersion = 0.1
fData(data)$dispFitted = 0.1
fData(data)$testable = ifelse((rowSums(counts(data) > 100)), TRUE, FALSE) #only test exon usage for genes with more than 100 counts

dataNorm = testForDEU(data)
dataNorm = estimatelog2FoldChanges(data)
result = DEUresultTable(data)
table(result$padjust != NA)

plot(result$meanBase, result[, "log2fold(negative/positive)"], log = "x",
     col = ifelse(res1$padjust < 0.1, "red", "black"), ylim = c(-4,4), main = "CD133 MvsA")

keep = rowSums(counts(dataNorm) > 10)