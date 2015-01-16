# Draw coverage plot according to http://davetang.org/muse/2013/09/07/creating-a-coverage-plot-in-r/
library(gplots)
library(RColorBrewer)

setwd("~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/copyNumber/")
list.files()

totCov = read.delim("141212_bedLike_average800.txt")

# Extract chr 7
# chr7 = totCov[totCov$chromosome_name %in% 7,c(4:12)]
# plot(chr7$GIC_011, type="l",
#      ylim = c(-4,4),
#      main = "Coverage plot", ylab = "Copy number",
#      xlab = "Chromosome 7", col = 'blue', lwd=2.5)
# lines(chr7$GIC_020, lwd=2.5, col = 'red')
# 
# # Try the same fpr chromosome 10
# chr10 = totCov[totCov$chromosome_name %in% 10,c(4:12)]
# plot(chr10$GIC_011, type="l",
#      ylim = c(-4,4),
#      main = "Coverage plot", ylab = "Copy number",
#      xlab = "Chromosome 10", col = 'blue', lwd=2.5)
# lines(chr10$GIC_020, lwd=2.5, col = 'red')

########### Test the heatmp method instead ###############
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

# Make heat map with copyNumber info
smallCov= rbind(totCov[c(1:200),],totCov[c(20000:20200),])
matCov = t(as.matrix(smallCov[,c(7:12)]))

# Retreive the indicies of the first occurance of the chromosome names for use in colsep
colSep = match(as.factor(c(1,17, 18)), smallCov$chromosome_name)

jpeg(filename="150116_smallHeatTest.jpeg", height=210, width=297,units="mm",
     res=300)
heatmap.2(matCov, cexRow=1.5, main="Copy number as inferred by RNA-seq",
          colsep=colSep, sepcolor="black", sepwidth=c(2,2),
          keysize=1, trace="none", key.title="CNV",
          col=myPalette, density.info="none", dendrogram="row", 
          labRow=row.names(matCov), xlab="Position",
          offsetRow=c(1,1), margins=c(2,7.5))
dev.off()

########### Full heatmap ###############
matCov = t(as.matrix(totCov[,c(7:12)]))



colSep = match(as.factor(c(1:22, "X", "Y")), totCov$chromosome_name)
jpeg(filename="150116_bigHeatTest.jpeg", height=210, width=297,units="mm",
     res=300)
# heatmap.2(matCov, cexRow=1.5, main="Copy number as inferred by RNA-seq",
#           colsep=colSep, sepcolor="black", sepwidth=c(2,2),
#           keysize=1, trace="none", key.title="CNV",
#           col=myPalette, density.info="none", dendrogram="row", 
#           labRow=row.names(matCov), xlab="Position",
#           offsetRow=c(1,1), margins=c(2,7.5))
dev.off()

quit("no")