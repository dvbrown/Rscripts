getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/wgcna/')
library(WGCNA)
source('/Users/d.brown6/Documents/Rscripts/140508_coexpressionFunctions.R')
options(stringsAsFactors=F)
list.files()

load('140415_justTheAgilentData.RData')
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/wgcna/manualCorrelation/')

######################################## CD133 coexpressed Genes ################################################
cd133 = correlateGeneWithGEM(dat, 'PROM1')
#write.table(cd133, './results/140520_cd133Coexpression.txt', row.names=T, sep='\t')
#write.table(abs(cd133), './results/140520_cd133CoexpressionAbs.txt', row.names=T, sep='\t')


plotCoexpression(cd133, 'CD133')

cd133genes = cd133[abs(cd133[,1]) > 2*sd(cd133[,1]) & cd133[,4] < 0.05,] # Use twice the standard deviation and significantly correlated
cd133Square = makeSquareCoexpressionMatrix(cd133genes, dat)

cd133Dissim = makeDissimilarity(cd133Square)

# Make a heatmap and store the colors of the idenitfied submodules
cd133Color = buildHeatMap(cd133Dissim, 'CD133')
# Make MDS plot
#pdf('cd133_MDS.pdf', paper='a4')
makeMDS(cd133Dissim, cd133Color, 'CD133')
#dev.off()

# Export to cytoscape
#cd133_cytoscape = cytoScapeInput(1-cd133Dissim, cd133Color, 'CD133')
#cd133_cytoscape1 = cd133_cytoscape[[2]]

######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')
# write.table(cd44, './results/140520_cd44Coexpression.txt', row.names=T, sep='\t')

plotCoexpression(cd44, 'CD44')

# Subset the dataframe with correlation values for those with high correlation and significance
#cd44genes = cd44[cd44[,2] > 0.1 & cd44[,4] < 0.05,]
cd44genes = cd44[abs(cd44[,1]) > 2*sd(cd44[,1]) & cd44[,4] < 0.05,] # Use twice the standard deviation and significantly correlated
cd44Square = makeSquareCoexpressionMatrix(cd44genes, dat)

cd44Dissim = makeDissimilarity(cd44Square)

# Make a heatmap and store the colors of the idenitfied submodules
cd44Color = buildHeatMap(cd44Dissim, 'CD44')
# Make MDS plot
#pdf('cd44_MDS.pdf', paper='a4')
makeMDS(cd44Dissim, cd44Color, 'CD44')
#dev.off()

cd44_cytoscape = cytoScapeInput(1-cd44Dissim, cd44Color, 'CD44')
