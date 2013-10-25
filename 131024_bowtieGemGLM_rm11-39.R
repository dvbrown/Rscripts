library(edgeR)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/')
source('~/Documents/Rscripts/annotateEnsembIDs.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')
files = list.files(pattern='*.txt')

colors = c('lightgreen','springgreen','cyan', 'lightblue')
dm = read.csv('designMatrix.csv')
dm = dm[c(2,3,4,6),]

f = lapply(files, read.delim, header=FALSE)
df1 = cbind(f[[2]],f[[3]],f[[4]],f[[6]])
df = df1[,c(2,4,6,8)]
row.names(df) = df1[,1]
colnames(df) = c('GIC_020', 'GIC_034', 'GIC_035', 'GIC_041')
labels = c('#020', '#034', '#035', '#041')

noFeatures = tail(df)
noCount = rownames(df) %in% c("no_feature","ambiguous","too_low_aQual",
                              "not_aligned","alignment_not_unique")
df = df[!noCount,]
totalCount = colSums(df)

#Build EdgeR objects
condition = c('long', 'long', 'short', 'short')
counts = DGEList(counts=df, group=condition)

# I remove genes with less than 0.5 cpm in 2 samples. For a library size of 20M this is 10 reads.
cpms = cpm(counts)
keep = rowSums(cpms >0.5) >=2
counts = counts[keep,]

#nomalise, plot MDS
d = calcNormFactors(counts)
plotMDS(d, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, main='MDS plot GIC RNA-seq batch1')
legend('topright', legend=c('Long-term','Short-term'), fill=c("darkgreen","blue"), cex=0.5)

boxplot(cpm(counts, log=T), main='Normalised counts RNA-seq batch1', ylab='Log2 CPM', col=colors, cex=1.25, las=2)

# Build the design matrix
design = model.matrix(~ libPrep + age + group, dm)

# Estimate dispersion. Cannot use GLM as lack degrees of freedom.
# Estimate using standrad empiral Bayes
d2 = estimateCommonDisp(d, design)

# Can only use common dispersion as lack replictes
#d2 = estimateTrendedDisp(d2, design)
#d2 = estimateTagwiseDisp(d2, design)

# Plot the dispersion. Tagwise vars is blue scatter. NB line is blue. Poisson line is black. Raw variance is maroon
plotMeanVar(d2, show.raw.vars=T, NBline=TRUE, main='Common dispersion GIC RNA-seq batch 1')
legend('topleft', legend=c('Poisson line', 'Neg Binomial line', 'Tagwise disp', 'Raw disp'), fill=c('black', 'steelblue', 'skyblue', 'maroon'), cex=0.8)

#plotBCV(d2, main='Biological variation GIC RNA-seq batch 1')

# Fit the GLM, interested in the group difference
f = glmFit(d2, design)
de = glmLRT(f, coef=4)

tt = topTags(de, n=nrow(d))
head(tt$table)

# Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot
rn = rownames(tt$table)
deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg, main='MA plot GIC RNA-seq batch1')

nc = cpm(d, normalized.lib.sizes=TRUE)
head(nc[rn,order(labels)],5)

result = ensembl_2_geneName(tt$table)
result = sort.dataframe(result, 8, highFirst=FALSE)
cutoff = result[result$FDR < 0.05,]
write.table(cpms,'GLMedgeR/remove011and039/131024_normalisedCPM',sep='\t')
write.table(result, './GLMedgeR/remove011and039/131024_shortVSlong.txt', sep='\t')
write.table(cutoff, './GLMedgeR/remove011and039/131024_shortVSlongDEgenes.txt', sep='\t')