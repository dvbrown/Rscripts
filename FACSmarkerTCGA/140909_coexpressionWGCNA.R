getwd()
library(WGCNA)
library(sqldf)
library(ggplot2)

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
######################################## CD44 coexpressed Genes ################################################
cd44 = correlateGeneWithGEM(dat, 'CD44')

# Use twice the standard deviation and significantly correlated
cd44genes = cd44[cd44[,1] > 2*sd(cd44[,1]) & cd44[,4] < 0.05,]

dbWriteTable(conn = db, name = "cd44Allgenes", value = as.data.frame(cd44), row.names = TRUE)
dbWriteTable(conn = db, name = "cd44CuttOff", value = as.data.frame(cd44genes), row.names = TRUE)
# cd44_cytoscape = cytoScapeInput(1-cd44Dissim, cd44Color, coexpressedShortList=cd44genes, 'CD44')

#### Make a histogram ####
cd133[,2] = "CD133"
cd44[,2] = "CD44"
df = rbind(cd133[,c(1,2)], cd44[,c(1,2)])
row.names(df) = NULL
df = as.data.frame(df)
colnames(df) = c("correlation", "Marker")
df = df[!is.na(df$correlation),]

df$correlation = as.character(df$correlation)
df$correlation = as.numeric(df$correlation)

ggplot(df, aes(x=correlation, fill=Marker)) + geom_density(alpha=.2) +
    xlab("Pearson correlation") + ylab("Frequency") + # Set axis labels
    ggtitle("Comparison of CD133 and CD44 correlation scores") +  # Set title
    theme_bw(base_size=18)

rm(cd44, cd44genes, cd133, cd133genes)
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

# Subset the dataframe with correlation values for those with high correlation and significance
L1CAMgenes = L1CAM[L1CAM[,1] > 2*sd(L1CAM[,1]) & L1CAM[,4] < 0.05,]

dbWriteTable(conn = db, name = "l1camAllgenes", value = as.data.frame(L1CAM), row.names = TRUE)
dbWriteTable(conn = db, name = "l1camCuttOff", value = as.data.frame(L1CAMgenes), row.names = TRUE)

# L1CAM_cytoscape = cytoScapeInput(1-L1CAMDissim, coexpressedShortList=L1CAMgenes, 'L1CAM')
rm(L1CAM, L1CAMgenes)

######################################## Olig2 coexpressed Genes ################################################
olig2 = correlateGeneWithGEM(dat, 'OLIG2')
olig2['OLIG2',]

# Subset the dataframe with correlation values for those with high correlation and significance
olig2genes = olig2[olig2[,1] > 2*sd(olig2[,1]) & olig2[,4] < 0.05,]

dbWriteTable(conn = db, name = "olig2Allgenes", value = as.data.frame(olig2), row.names = TRUE)
dbWriteTable(conn = db, name = "olig2CuttOff", value = as.data.frame(olig2genes), row.names = TRUE)

rm(olig2, olig2genes)
######################################## GFAP coexpressed Genes ################################################
gfap = correlateGeneWithGEM(dat, 'GFAP')
gfap['GFAP',]

# Subset the dataframe with correlation values for those with high correlation and significance
gfapgenes = gfap[gfap[,1] > 2*sd(gfap[,1]) & gfap[,4] < 0.05,]

dbWriteTable(conn = db, name = "gfapAllgenes", value = as.data.frame(gfap), row.names = TRUE)
dbWriteTable(conn = db, name = "gfapCuttOff", value = as.data.frame(gfapgenes), row.names = TRUE)

rm(gfap, gfapgenes)
######################################## YKL40 coexpressed Genes ################################################
ykl40 = correlateGeneWithGEM(dat, 'CHI3L1')
ykl40['CHI3L1',]

# Subset the dataframe with correlation values for those with high correlation and significance
ykl40genes = ykl40[ykl40[,1] > 2*sd(ykl40[,1]) & ykl40[,4] < 0.05,]

dbWriteTable(conn = db, name = "ykl40Allgenes", value = as.data.frame(ykl40), row.names = TRUE)
dbWriteTable(conn = db, name = "ykl40CuttOff", value = as.data.frame(ykl40genes), row.names = TRUE)

rm(ykl40, ykl40genes)
######################################## SOX2 coexpressed Genes ################################################
sox2 = correlateGeneWithGEM(dat, 'SOX2')
sox2['SOX2',]
# Subset the dataframe with correlation values for those with high correlation and significance
sox2genes = sox2[sox2[,1] > 2*sd(sox2[,1]) & sox2[,4] < 0.05,]

dbWriteTable(conn = db, name = "sox2Allgenes", value = as.data.frame(sox2), row.names = TRUE)
dbWriteTable(conn = db, name = "sox2CuttOff", value = as.data.frame(sox2genes), row.names = TRUE)

rm(sox2, sox2genes)

######################################## B-3- tubulin coexpressed Genes ################################################
b3tub = correlateGeneWithGEM(dat, 'TUBB3')
b3tub['TUBB3',]
# Subset the dataframe with correlation values for those with high correlation and significance
b3tubgenes = b3tub[b3tub[,1] > 2*sd(b3tub[,1]) & b3tub[,4] < 0.05,]

dbWriteTable(conn = db, name = "tubb3Allgenes", value = as.data.frame(b3tub), row.names = TRUE)
dbWriteTable(conn = db, name = "tubb3CuttOff", value = as.data.frame(b3tubgenes), row.names = TRUE)

rm(b3tub, b3tubgenes)

dbDisconnect(db)