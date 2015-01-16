# Draw coverage plot according to http://davetang.org/muse/2013/09/07/creating-a-coverage-plot-in-r/

setwd("~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/copyNumber/")
list.files()

totCov = read.delim("141212_bedLike_average800.txt")

# Extract chr 7
smallCov = totCov[totCov$chromosome_name %in% 7,c(4:12)]
plot(smallCov$GIC_011, type="l",
     ylim = c(-4,4),
     main = "Coverage plot", ylab = "Copy number",
     xlab = "Chromosome 7", col = 'blue', lwd=2.5)
lines(smallCov$GIC_020, lwd=2.5, col = 'red')

# Try the same fpr chromosome 10
chr10 = totCov[totCov$chromosome_name %in% 10,c(4:12)]
plot(chr10$GIC_011, type="l",
     ylim = c(-4,4),
     main = "Coverage plot", ylab = "Copy number",
     xlab = "Chromosome 10", col = 'blue', lwd=2.5)
lines(chr10$GIC_020, lwd=2.5, col = 'red')

########### Try out the heatmp method instead ###############