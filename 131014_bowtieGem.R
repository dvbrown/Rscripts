library(edgeR)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/')
files = list.files(pattern='*.txt')

f = lapply(files, read.delim, header=FALSE)
df1 = cbind(f[[1]],f[[2]],f[[3]],f[[4]],f[[5]],f[[6]])
df = df1[,c(2,4,6,8,10,12)]
row.names(df) = df1[,1]
colnames(df) = c('GIC_011', 'GIC_020', 'GIC_034', 'GIC_035', 'GIC_039', 'GIC_041')
labels = c('#011', '#020', '#034', '#035', '#039', '#041')

noFeatures = tail(df)
df = df[c(1:203288),]
totalCount = colSums(df)

#Build EdgeR objects
condition = c('long', 'long', 'long', 'short', 'short', 'short')
counts = DGEList(counts=df, group=condition)

#critical step In edgeR, it is recommended to remove features without at least 1 read per million in n of the samples,
#where n is the size of the smallest group of replicates (here, n = 3 for the knockdown group)
cpms = cpm(counts)
keep = rowSums(cpms >1) >=3
counts = counts[keep,]

#nomalise, plot MDS
d = calcNormFactors(counts)
plotMDS(d, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, main='MDS plot GIC RNA-seq batch1')

d = estimateCommonDisp(d)
#d = estimateTagwiseDisp(d)

# Plot the dispersions
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
plotBCV(d)

# Differential expression testing
de = exactTest(d, pair=c("long","short"))
tt = topTags(de, n=nrow(d))
head(tt$table)

# Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot
deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg)
