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

# Use twice the standard deviation and significantly correlated
cd133genes = cd133[(cd133[,1]) > 2*sd(cd133[,1]) & cd133[,4] < 0.05,]
# write.table(cd133genes, './140527_cd133Cutoff.txt', sep='\t')
cd133Square = makeSquareCoexpressionMatrix(cd133genes, dat)

cd133Dissim = makeDissimilarity(cd133Square)

# Make a heatmap
pdf(file="cd133_corr_HeatMap.pdf", paper="a4")
buildCorrelationHeatMap(dat, cd133Dissim, 'CD133')
dev.off()

pdf(file="cd133_tom_HeatMap.pdf", paper="a4")
buildTOMHeatMap(dat, cd133Dissim, "CD133")
dev.off()

# Make MDS plot
#pdf('cd133_MDS.pdf', paper='a4')
moduleCol = flashClust(as.dist(cd133Dissim))
makeMDS(cd133Dissim, 'CD133')
#dev.off()

# Export to cytoscape
cd133_cytoscape = cytoScapeInput(1-cd133Dissim,coexpressedShortList=cd133genes, 'CD133')

rm(cd133, cd133genes)
######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')
# write.table(cd44, './results/140526_cd44Coexpression.txt', row.names=T, sep='\t')

plotCoexpression(cd44, 'CD44')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
cd44genes = cd44[cd44[,1] > 3*sd(cd44[,1]) & cd44[,4] < 0.05,]
write.table(cd44genes, './140529_cd44Cutoff.txt', sep='\t')

cd44Square = makeSquareCoexpressionMatrix(cd44genes, dat)
cd44Dissim = makeDissimilarity(cd44Square)

# Make a heatmap
pdf(file="cd44_corr_HeatMap.pdf", paper="a4")
buildCorrelationHeatMap(dat, cd44Dissim, 'CD44')
dev.off()

pdf(file="cd44_tom_HeatMap.pdf", paper="a4")
buildTOMHeatMap(dat, cd44Dissim, "CD44")
dev.off()

makeMDS(cd44Dissim, 'CD44')
#dev.off()

cd44_cytoscape = cytoScapeInput(1-cd44Dissim, cd44Color, coexpressedShortList=cd44genes, 'CD44')

rm(cd44, cd44genes)
######################################## CD15 coexpressed Genes ################################################
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/CD15')
cd15 = correlateGeneWithGEM(dat, 'FUT4')
write.table(cd15, './140529_cd15Coexpression.txt', row.names=T, sep='\t')

plotCoexpression(cd15, 'CD15')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
cd15genes = cd15[cd15[,1] > 3*sd(cd15[,1]) & cd15[,4] < 0.05,]
write.table(cd15genes, './140529_cd15Cutoff.txt', sep='\t')

cd15Square = makeSquareCoexpressionMatrix(cd15genes, dat)
cd15Dissim = makeDissimilarity(cd15Square)

# Make a heatmap
pdf(file="cd15_corr_HeatMap.pdf", paper="a4")
buildCorrelationHeatMap(dat, cd15Dissim, 'CD15')
dev.off()

pdf(file="cd15_tom_HeatMap.pdf", paper="a4")
buildTOMHeatMap(dat, cd15Dissim, "CD15")
dev.off()

#pdf('cd44_MDS.pdf', paper='a4')
makeMDS(cd15Dissim, 'CD15')
#dev.off()

cd15_cytoscape = cytoScapeInput(1-cd15Dissim, coexpressedShortList=cd15genes, 'CD15')

# Remove CD related objects from the environment
rm(cd15, cd15genes)
######################################## Integrin alpha 6 coexpressed Genes ################################################
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ITGA6/')

ITGA6 = correlateGeneWithGEM(dat, 'ITGA6')
write.table(ITGA6, './140529_ITGA6Coexpression.txt', row.names=T, sep='\t')

plotCoexpression(ITGA6, 'ITGA6')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
ITGA6genes = ITGA6[ITGA6[,1] > 3*sd(ITGA6[,1]) & ITGA6[,4] < 0.05,]
write.table(ITGA6genes, './140529_ITGA6Cutoff.txt', sep='\t')

ITGA6Square = makeSquareCoexpressionMatrix(ITGA6genes, dat)
ITGA6Dissim = makeDissimilarity(ITGA6Square)

# Make a heatmap
pdf(file="ITGA6_corr_HeatMap.pdf", paper="a4")
buildCorrelationHeatMap(dat, ITGA6Dissim, 'ITGA6')
dev.off()

pdf(file="ITGA6_tom_HeatMap.pdf", paper="a4")
buildTOMHeatMap(dat, ITGA6Dissim, "ITGA6")
dev.off()

#pdf('cd44_MDS.pdf', paper='a4')
makeMDS(ITGA6Dissim, 'ITGA6')
#dev.off()

ITGA6_cytoscape = cytoScapeInput(1-ITGA6Dissim, coexpressedShortList=ITGA6genes, 'ITGA6')

rm(ITGA6, ITGA6genes)
######################################## Aldehyde dehydorgenase 1 coexpressed Genes ################################################
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/ALDH1/')

ALDH1 = correlateGeneWithGEM(dat, 'ALDH1A1')
write.table(ALDH1, './140529_ALDH1Coexpression.txt', row.names=T, sep='\t')

plotCoexpression(ALDH1, 'ALDH1')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
ALDH1genes = ALDH1[ALDH1[,1] > 3*sd(ALDH1[,1]) & ALDH1[,4] < 0.05,]
write.table(ALDH1genes, './140529_ALDH1Cutoff.txt', sep='\t')

ALDH1Square = makeSquareCoexpressionMatrix(ALDH1genes, dat)
ALDH1Dissim = makeDissimilarity(ALDH1Square)

# Make a heatmap
pdf(file="ALDH1_corr_HeatMap.pdf", paper="a4")
buildCorrelationHeatMap(dat, ALDH1Dissim, 'ALDH1')
dev.off()

pdf(file="ALDH1_tom_HeatMap.pdf", paper="a4")
buildTOMHeatMap(dat, ALDH1Dissim, "ALDH1")
dev.off()

pdf('ALDH1_MDS.pdf', paper='a4')
makeMDS(ALDH1Dissim, 'ALDH1')
dev.off()

ALDH1_cytoscape = cytoScapeInput(1-ALDH1Dissim, coexpressedShortList=ALDH1genes, 'ALDH1')

rm(ALDH1, ALDH1genes)
######################################## L1CAM coexpressed Genes ################################################
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/L1CAM/')

L1CAM = correlateGeneWithGEM(dat, 'L1CAM')
write.table(L1CAM, './140529_L1CAMCoexpression.txt', row.names=T, sep='\t')

plotCoexpression(L1CAM, 'L1CAM')

# Subset the dataframe with correlation values for those with high correlation and significance
# Use twice the standard deviation and significantly correlated
L1CAMgenes = L1CAM[L1CAM[,1] > 3*sd(L1CAM[,1]) & L1CAM[,4] < 0.05,]
write.table(L1CAMgenes, './140529_L1CAMCutoff.txt', sep='\t')

L1CAMSquare = makeSquareCoexpressionMatrix(L1CAMgenes, dat)
L1CAMDissim = makeDissimilarity(L1CAMSquare)

# Make a heatmap
pdf(file="L1CAM_corr_HeatMap.pdf", paper="a4")
buildCorrelationHeatMap(dat, L1CAMDissim, 'L1CAM')
dev.off()

pdf(file="L1CAM_tom_HeatMap.pdf", paper="a4")
buildTOMHeatMap(dat, L1CAMDissim, "L1CAM")
dev.off()

pdf('L1CAM_MDS.pdf', paper='a4')
makeMDS(L1CAMDissim, 'L1CAM')
dev.off()

L1CAM_cytoscape = cytoScapeInput(1-L1CAMDissim, coexpressedShortList=L1CAMgenes, 'L1CAM')
rm(L1CAM, L1CAMgenes)
################################################################
#doublePos = corAndPvalue(x=dat[,c("PROM1", "CD44")], y=dat)