getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/wgcna/')
library(WGCNA)
source('/Users/d.brown6/Documents/Rscripts/140508_coexpressionFunctions.R')
options(stringsAsFactors=F)
list.files()

load('140415_justTheAgilentData.RData')

######################################## CD133 coexpressed Genes ################################################
cd133 = correlateGeneWithGEM(dat, 'PROM1')

plotCoexpression(cd133, 'CD133')

cd133genes = cd133[abs(cd133[,1]) > 2*sd(cd133[,1]) & cd133[,4] < 0.05,] # Use twice the standard deviation and significantly correlated





######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')

plotCoexpression(cd44, 'CD44')

# Subset the dataframe with correlation values for those with high correlation and significance
#cd44genes = cd44[cd44[,2] > 0.1 & cd44[,4] < 0.05,]
cd44genes = cd44[abs(cd44[,1]) > 2*sd(cd44[,1]) & cd44[,4] < 0.05,] # Use twice the standard deviation and significantly correlated