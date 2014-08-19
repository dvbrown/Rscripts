getwd()
library(WGCNA)
library(ggplot2)
source('/Users/d.brown6/Documents/Rscripts/140508_coexpressionFunctions.R')
options(stringsAsFactors=F)
list.files()

dat = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')

######################################## CD133 coexpressed Genes ################################################
cd133 = correlateGeneWithGEM(dat, 'PROM1')
# write.table(cd133, './140526_cd133Coexpression.txt', row.names=T, sep='\t')
plotCoexpression(cd133, 'CD133')

# Use twice the standard deviation and significantly correlated
cd133genes = cd133[(cd133[,1]) > 2*sd(cd133[,1]) & cd133[,4] < 0.05,]
# write.table(cd133genes, './140527_cd133Cutoff.txt', sep='\t')
cd133Square = makeSquareCoexpressionMatrix(cd133genes, dat)

cd133Dissim = makeDissimilarity(cd133Square)

######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')
plotCoexpression(cd44, 'CD44')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
cd44genes = cd44[cd44[,1] > 3*sd(cd44[,1]) & cd44[,4] < 0.05,]

cd44Square = makeSquareCoexpressionMatrix(cd44genes, dat)
cd44Dissim = makeDissimilarity(cd44Square)

######################################## CD15 coexpressed Genes ################################################
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15')
cd15 = correlateGeneWithGEM(dat, 'FUT4')

plotCoexpression(cd15, 'CD15')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
cd15genes = cd15[cd15[,1] > 3*sd(cd15[,1]) & cd15[,4] < 0.05,]
