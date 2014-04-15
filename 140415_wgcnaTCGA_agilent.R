getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
library(WGCNA)
library(biomaRt)
options(stringsAsFactors=F)
list.files()

#################################### Data input and clustering visualisations ##############################################
data = read.delim('140110_agilentNoNulls.txt')

# We now transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(data))
dat = as.matrix(datExpr0)
save.image('140415_justTheAgilentData.RData')
########################################################## Standard correlations ###############################################
load('140415_justTheAgilentData.RData')

# Calculate the correlation between PROM1 expression and all the genes in TCGA GBM
prom1CorrPval = corAndPvalue(x=datExpr0[,'PROM1'], y=dat)

#Extract the correlation and p-value from the returned list
prom1C = prom1CorrPval$cor
prom1Cpower = prom1C^2

prom1P = prom1CorrPval$p
prom1FDR = p.adjust(prom1P, method='fdr')

# par(mfrow=c(2,1))
# hist(prom1Cpower, main='Prom1 correlations', breaks='FD', xlab='Weighted correlation values')
# hist(prom1FDR, main='Prom1 p-values', breaks='FD', xlab='FDR corrected p-values')
# par(mfrow=c(1,1))

result = t(rbind(prom1C, prom1Cpower, prom1P, prom1FDR))
colnames(result) = c('correlation', 'weighted_correlation', 'p-value', 'FDR')
write.table(result, './wgcna/140415_standardCorrelation.txt', sep='\t')

par(mfrow=c(2,2))
hist(prom1C, main='Prom1 correlations', breaks='FD', xlab='Raw correlation values')
hist(prom1Cpower, main='Prom1 weighted correlations', breaks='FD', xlab='R squared values')
hist(prom1P, main='Prom1 p-values', breaks='FD', xlab='Raw p-values')
hist(prom1FDR, main='Prom1FDR corrected p-values', breaks='FD', xlab='FDR corrected p-values')
par(mfrow=c(1,1))

qqnorm(result[,1], main='Distrubution of PROM1 correlated genes')
qqline(result[,1])

################################################################################################################### 


##################################### Build a simlarity matrix ######################################################