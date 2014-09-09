getwd()
library(WGCNA)
library(sqldf)

dat = read.delim("~/Documents/public-datasets/cancerBrowser/deDupAgilent/140526_agilentDedupPatients.txt", row.names=1)
setwd('/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/results/')
source('~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')
list.files()

#########################################  Initalise a connection to sqlite database ###################################################################
db <- dbConnect(SQLite(), dbname="/Users/d.brown6/Documents/public-datasets/cancerBrowser/deDupAgilent/coexpression.sqlite")

######################################## CD133 coexpressed Genes ################################################
cd133 = correlateGeneWithGEM(dat, 'PROM1')

# Use twice the standard deviation and significantly correlated
cd133genes = cd133[(cd133[,1]) > 2*sd(cd133[,1]) & cd133[,4] < 0.05,]

dbWriteTable(conn = db, name = "cd133Allgenes", value = as.data.frame(cd133), row.names = TRUE)
dbWriteTable(conn = db, name = "cd133CuttOff", value = as.data.frame(cd133genes), row.names = TRUE)

# Export to cytoscape
# cd133_cytoscape = cytoScapeInput(1-cd133Dissim,coexpressedShortList=cd133genes, 'CD133')
rm(cd133, cd133genes)
######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')

# Use twice the standard deviation and significantly correlated
cd44genes = cd44[cd44[,1] > 2*sd(cd44[,1]) & cd44[,4] < 0.05,]

dbWriteTable(conn = db, name = "cd44Allgenes", value = as.data.frame(cd44), row.names = TRUE)
dbWriteTable(conn = db, name = "cd44CuttOff", value = as.data.frame(cd44genes), row.names = TRUE)
# cd44_cytoscape = cytoScapeInput(1-cd44Dissim, cd44Color, coexpressedShortList=cd44genes, 'CD44')

rm(cd44, cd44genes)
######################################## CD15 coexpressed Genes ################################################
cd15 = correlateGeneWithGEM(dat, 'FUT4')

# Use twice the standard deviation and significantly correlated
cd15genes = cd15[cd15[,1] > 2*sd(cd15[,1]) & cd15[,4] < 0.05,]

dbWriteTable(conn = db, name = "cd15Allgenes", value = as.data.frame(cd15), row.names = TRUE)
dbWriteTable(conn = db, name = "cd15CuttOff", value = as.data.frame(cd15genes), row.names = TRUE)
#cd15_cytoscape = cytoScapeInput(1-cd15Dissim, coexpressedShortList=cd15genes, 'CD15')

# Remove CD related objects from the environment
rm(cd15, cd15genes)
######################################## Integrin alpha 6 coexpressed Genes ################################################
ITGA6 = correlateGeneWithGEM(dat, 'ITGA6')

# Use twice the standard deviation and significantly correlated
ITGA6genes = ITGA6[ITGA6[,1] > 2*sd(ITGA6[,1]) & ITGA6[,4] < 0.05,]

dbWriteTable(conn = db, name = "itga6Allgenes", value = as.data.frame(ITGA6), row.names = TRUE)
dbWriteTable(conn = db, name = "itag6CuttOff", value = as.data.frame(ITGA6genes), row.names = TRUE)
# ITGA6_cytoscape = cytoScapeInput(1-ITGA6Dissim, coexpressedShortList=ITGA6genes, 'ITGA6')

rm(ITGA6, ITGA6genes)
######################################## Aldehyde dehydorgenase 1 coexpressed Genes ################################################
ALDH1 = correlateGeneWithGEM(dat, 'ALDH1A1')

# Subset the dataframe with correlation values for those with high correlation and significance
ALDH1genes = ALDH1[ALDH1[,1] > 2*sd(ALDH1[,1]) & ALDH1[,4] < 0.05,]

dbWriteTable(conn = db, name = "aldh1Allgenes", value = as.data.frame(ALDH1), row.names = TRUE)
dbWriteTable(conn = db, name = "aldh1CuttOff", value = as.data.frame(ALDH1genes), row.names = TRUE)
# ALDH1_cytoscape = cytoScapeInput(1-ALDH1Dissim, coexpressedShortList=ALDH1genes, 'ALDH1')

rm(ALDH1, ALDH1genes)
######################################## L1CAM coexpressed Genes ################################################
L1CAM = correlateGeneWithGEM(dat, 'L1CAM')

plotCoexpression(L1CAM, 'L1CAM')

# Subset the dataframe with correlation values for those with high correlation and significance
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