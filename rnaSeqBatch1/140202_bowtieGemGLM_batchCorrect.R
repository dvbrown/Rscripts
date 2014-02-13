library(edgeR)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/')
source('~/Documents/Rscripts/annotateEnsembIDs.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')
files = list.files(pattern='*.txt')

colors = c('darkgreen', 'lightgreen','springgreen','cyan', 'blue1', 'lightblue')
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
cpmLog = cpm(counts, log=T, prior.count=0.5)
keep = rowSums(cpms >0.5) >=3
counts = counts[keep,]

#nomalise, plot MDS
d = calcNormFactors(counts)

############################################## MDS plots and removing batch effects #########################################
par(mfrow=c(2,2))
# Plots according to the biological coefficient of variation
plotMDS(d, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, 
       main='MDS plot GIC RNA-seq batch1 original')
legend('topright', legend=c('Long-term','Short-term'), fill=c("darkgreen","blue"), cex=0.33)
# 
facs = removeBatchEffect(cpms, batch=dm$facsSort)
# Retreive cpms from FACS
logCpmFacs  = removeBatchEffect(cpmLog, batch=dm$facsSort)

plotMDS(facs, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, 
        main='MDS plot GIC RNA-seq batch1 facs sorts')

lib = removeBatchEffect(d, batch=dm$libPrep)
plotMDS(lib, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, 
        main='MDS plot GIC RNA-seq batch1 library prep')

doubleCorrect = removeBatchEffect(d, batch=dm$facsSort, batch2=dm$libPrep)
plotMDS(doubleCorrect, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, 
        main='MDS plot GIC RNA-seq batch1 \nlibrary and FACS correct')

# plotMDS(d, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, 
#         main='MDS plot GIC RNA-seq batch1', dim.plot=c(2,3))
# #legend('topright', legend=c('Long-term','Short-term'), fill=c("darkgreen","blue"), cex=0.33)
# 
# plotMDS(d, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, 
#         main='MDS plot GIC RNA-seq batch1', dim.plot=c(3,4))
# #legend('topright', legend=c('Long-term','Short-term'), fill=c("darkgreen","blue"), cex=0.33)
# 
# plotMDS(d, labels=labels, col = c("darkgreen","blue")[factor(condition)], cex=1.25, 
#         main='MDS plot GIC RNA-seq batch1', dim.plot=c(4,5))
#legend('topright', legend=c('Long-term','Short-term'), fill=c("darkgreen","blue"), cex=0.33)

##############################################################################################################################

boxplot(cpm(counts, log=T), main='Normalised counts RNA-seq batch1', ylab='Log2 CPM', col=colors, cex=1.25, las=2)

# Build the design matrix
design = model.matrix(~ facsSort + age + group, dm)

# Estimate dispersion, adjusting for the GLM
d2 = estimateGLMCommonDisp(d, design, verbose=T)
# Condense the disersion estimate by estiamting a trend for count size
d2 = estimateGLMTrendedDisp(d2, design)
# Empriical Bayes tagwise dispersion
d2 = estimateGLMTagwiseDisp(d2, design)

# Plot the dispersions. Tagwise vars is blue scatter. NB line is blue. Poisson line is black. Raw variance is maroon
plotMeanVar(d2, show.tagwise.vars=TRUE, NBline=TRUE, main='Fitted dispersion GIC RNA-seq batch 1')
legend('topleft', legend=c('Poisson line', 'Neg Binomial line', 'Tagwise disp', 'Raw disp'), fill=c('black', 'steelblue', 'skyblue', 'maroon'), cex=0.8)


plotBCV(d2, main='Biological variation GIC RNA-seq batch 1')

# Fit the GLM, interested in the group difference
# Fit the GLM per gene
f = glmFit(d2, design)
de = glmLRT(f, coef=4)

tt = topTags(de, n=nrow(d))
head(tt$table)
summary(decideTestsDGE(de))

# The quasi likihood metod doesn't give any differentially expressed genes
# fit <- glmQLFTest(d2,design,robust=TRUE,plot=TRUE)
# t = topTags(fit, n=nrow(d))
# head(t$table)
# summary(decideTestsDGE(fit))

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
write.table(cpmFacs,'GLMedgeR/140203_facsBatch/140213_normalisedCPM_facs.txt',sep='\t')
write.table(logCpmFacs,'GLMedgeR/140203_facsBatch/140213_normalisedLog_CPM_facs.txt',sep='\t')
write.table(result, './GLMedgeR/140203_facsBatch/140203_shortVSlong.txt', sep='\t')
write.table(cutoff, './GLMedgeR/140203_facsBatch/140203_shortVSlongDEgenes.txt', sep='\t')
write.table(cutoffLib, './GLMedgeR/140203_facsBatch/140203_shortVSlongLiberalDE.txt', sep='\t')
