library(edgeR)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/')
source('~/Documents/Rscripts/annotateEnsembIDs.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')


colors = c('darkgreen', 'lightgreen','springgreen','cyan', 'blue1', 'lightblue')
dm = read.csv('designMatrix.csv')

#### Modify the design matrix to include other groupings from the PCA and hirechaical clustering ####
pca = c("left", "left", "right", "right", "right", "left")
hire = c("right", "left", "right", "right", "left", "left")

#dm$group = pca
dm$group = hire
####

df = read.delim("131021_normalisedCPM.txt")
colnames(df) = c('MU011', 'MU020', 'MU034', 'MU035', 'MU039', 'GIC_MU041')
labels = c('MU011', 'MU020', 'MU034', 'MU035', 'MU039', 'GIC_MU041')

noFeatures = tail(df)
noCount = rownames(df) %in% c("no_feature","ambiguous","too_low_aQual",
                              "not_aligned","alignment_not_unique")
df = df[!noCount,]
totalCount = colSums(df)

#Build EdgeR objects
# condition = c('long', 'long', 'long', 'short', 'short', 'short')
counts = DGEList(counts=df, group=dm$group)

# I remove genes with less than 0.5 cpm in 3 samples. For a library size of 20M this is 10 reads.
cpms = cpm(counts)
cpmLog = cpm(counts, log=T, prior.count=0.5)
keep = rowSums(cpms >0.5) >=3
counts = counts[keep,]

#nomalise, plot MDS
d = calcNormFactors(counts)

# Build the design matrix
design = model.matrix(~ facsSort + age + group, dm)

# Estimate dispersion, adjusting for the GLM
d2 = estimateGLMCommonDisp(d, design, verbose=T)
# Condense the disersion estimate by estiamting a trend for count size
d2 = estimateGLMTrendedDisp(d2, design)
# Empriical Bayes tagwise dispersion
d2 = estimateGLMTagwiseDisp(d2, design)

# Fit the GLM, interested in the group difference
# Fit the GLM per gene
f = glmFit(d2, design)
de = glmLRT(f, coef=4)

tt = topTags(de, n=nrow(d))
head(tt$table)
summary(decideTestsDGE(de))

# Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot
rn = rownames(tt$table)
deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg, main='MA plot GIC RNA-seq batch1')
abline(h=c(-1, 1), col="blue")

nc = cpm(d, normalized.lib.sizes=TRUE)
head(nc[rn,order(labels)],5)

# Have a look at the counts of the most significant DE genes. Look for consistency in the counts
o = order(de$table$PValue)
cpm(d2)[o[1:50],]

result = ensembl_2_geneName(tt$table)
result = sort.dataframe(result, 8, highFirst=FALSE)
cutoff = result[result$FDR < 0.05,]
cutoffLib = result[result$FDR < 0.1 & abs(result$logFC) > 1,]

########################################################### Write out results ############################################## 
#write.table(cpmFacs,'GLMedgeR/140203_facsBatch/140213_normalisedCPM_facs.txt',sep='\t')
#write.table(logCpmFacs,'GLMedgeR/140203_facsBatch/140213_normalisedLog_CPM_facs.txt',sep='\t')
#write.table(result, './GLMedgeR/140203_facsBatch/140203_shortVSlong.txt', sep='\t')
#write.table(cutoff, './GLMedgeR/140203_facsBatch/140203_shortVSlongDEgenes.txt', sep='\t')
#write.table(cutoffLib, './GLMedgeR/140203_facsBatch/140203_shortVSlongLiberalDE.txt', sep='\t')
