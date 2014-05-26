getwd()
library(WGCNA)
source('/Users/d.brown6/Documents/Rscripts/140508_coexpressionFunctions.R')
options(stringsAsFactors=F)
list.files()

dat = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')

######################################## CD133 coexpressed Genes ################################################
cd133 = correlateGeneWithGEM(dat, 'PROM1')
# write.table(cd133, './140526_cd133Coexpression.txt', row.names=T, sep='\t')
plotCoexpression(cd133, 'CD133')

# subsample and return summary statistics
# cd133SubsamplesCorr = subsample10times(dat, "PROM1", 100, "correlation")
# cd133SubsamplesFDR = subsample10times(dat, "PROM1", 100, "FDR")
# plotResampling(cd133SubsamplesCorr, cd133SubsamplesFDR, cd133, "PROM1")
# rm(cd133SubsamplesCorr, cd133SubsamplesFDR)

cd133genes = cd133[(cd133[,1]) > 3*sd(cd133[,1]) & cd133[,4] < 0.05,] # Use twice the standard deviation and significantly correlated
cd133Square = makeSquareCoexpressionMatrix(cd133genes, dat)

cd133Dissim = makeDissimilarity(cd133Square)

# Make a heatmap and store the colors of the idenitfied submodules
cd133Color = buildHeatMap(cd133Dissim, 'CD133')
# Make MDS plot
#pdf('cd133_MDS.pdf', paper='a4')
makeMDS(cd133Dissim, cd133Color, 'CD133')
#dev.off()

# Export to cytoscape
cd133_cytoscape = cytoScapeInput(1-cd133Dissim, cd133Color,coexpressedShortList=cd133genes, 'CD133')

######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')
# write.table(cd44, './results/140526_cd44Coexpression.txt', row.names=T, sep='\t')

plotCoexpression(cd44, 'CD44')

# Subset the dataframe with correlation values for those with high correlation and significance
#cd44genes = cd44[cd44[,2] > 0.1 & cd44[,4] < 0.05,]
cd44genes = cd44[(cd44[,1]) > 3*sd(cd44[,1]) & cd44[,4] < 0.05,] # Use twice the standard deviation and significantly correlated
cd44Square = makeSquareCoexpressionMatrix(cd44genes, dat)

cd44Dissim = makeDissimilarity(cd44Square)

# Make a heatmap and store the colors of the idenitfied submodules
cd44Color = buildHeatMap(cd44Dissim, 'CD44')
# Make MDS plot
#pdf('cd44_MDS.pdf', paper='a4')
makeMDS(cd44Dissim, cd44Color, 'CD44')
#dev.off()

cd44_cytoscape = cytoScapeInput(1-cd44Dissim, cd44Color, coexpressedShortList=cd44genes, 'CD44')

######################################## CD133 CD44 double positive ################################################
doublePos = intersect(row.names(cd133genes), row.names(cd44genes))

doublePos = corAndPvalue(x=dat[,c("PROM1", "CD44")], y=dat)