# Use 140508_coexpression as a guide to heatmappin

getwd()
library(WGCNA)
source('/Users/d.brown6/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')
options(stringsAsFactors=F)
list.files()

agilentGem = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
rnaseqGem = read.delim("~/Documents/public-datasets/cancerBrowser/TCGA_GBM_exp_HiSeqV2-2014-05-02/genomicMatrix", row.names=1)

setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')

######################################## CD133 coexpressed Genes, for Agilent ################################################
cd133 = correlateGeneWithGEM(dat, 'PROM1')

# Use twice the standard deviation and significantly correlated
cd133genes = cd133[(cd133[,1]) > 2*sd(cd133[,1]) & cd133[,4] < 0.05,]
# write.table(cd133genes, './140527_cd133Cutoff.txt', sep='\t')
cd133Square = makeSquareCoexpressionMatrix(cd133genes, dat)
