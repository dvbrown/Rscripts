library(edgeR)
setwd('~/Documents/RNAdata/danBatch1/batch1FeatureCounts/')
files = list.files(pattern='*.txt')

# Read in files as a list of dataframes
f = lapply(files, read.delim, skip=1, header=TRUE)
g011 = f[[1]][,c(1,7)]
g020 = f[[2]][,c(1,7)]
g034 = f[[3]][,c(1,7)]
g035 = f[[4]][,c(1,7)]
g039 = f[[5]][,c(1,7)]
g041 = f[[6]][,c(1,7)]

# Unlist into 1 dataframe
df = data.frame(g011[,2], g020[,2], g034[,2], g035[,2], g039[,2], g041[,2], row.names=g011[,1])
colnames(df) = c('GIC_011', 'GIC_020', 'GIC_034', 'GIC_035', 'GIC_039', 'GIC_041')
labels = c('#011', '#020', '#034', '#035', '#039', '#041')

# Total reads
totalCount = colSums(df)

#Build EdgeR objects
condition = c('long', 'long', 'long', 'short', 'short', 'short')
counts = DGEList(counts=df, group=condition)

# Critical step In edgeR, it is recommended to remove features without at least 1 read per million in n of the samples,
# Where n is the size of the smallest group of replicates (here, n = 1 for patients)
cpms = cpm(counts)
keep = rowSums(cpms >=1) >=1
counts = counts[keep,]

#nomalise, plot MDS
d = calcNormFactors(counts)
plotMDS(d, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, main='MDS plot GIC RNA-seq batch1')

#d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)

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