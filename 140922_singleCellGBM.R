library(sqldf)
library(GSVA)
library(gplots)

setwd('~/Documents/public-datasets/RNA-seq/anoop2014_singleCellGBM/')
db <- dbConnect(SQLite(), dbname='~/Documents/public-datasets/cancerBrowser/deDupAgilent/coexpression.sqlite')

subsetSamples <- function (dataFrame, sampleNameStub) {
    # The dataframe with the measurements
    # A string representing the basename of the sample eg MGH26
    sample = paste(sampleNameStub, '*', sep='_')
    samples = grep(sample, colnames(dataFrame), value=T)
    result = dataFrame[,samples]
    return (result)
}

measureSignatures <- function (dataFrame, signatureList) {
    # First argument is the dataframe containing the gene expression measurments
    # Second argument is the list containing the names of genes in the signature
    dataMatrix = as.matrix(dataFrame)
    sigScore = gsva(dataMatrix, signatureList,  rnaseq=F, verbose=T, parallel.sz=1)
    result = t(sigScore$es.obs)
    return (result)
}
############################################## IO #############################################

dbListTables(db)
cd133Sig = dbReadTable(db, "cd133CuttOff")
cd44Sig = dbReadTable(db, "cd44CuttOff")
cd15 = dbReadTable(db, "cd15CuttOff")
signatures = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig), 'CD15' = row.names(cd15))
rm(cd133Sig, cd44Sig, cd15)

data = read.delim('GSE57872_GBM_data_matrix.txt', row.names=1)
annotation = read.delim('sample.txt')

############################################### Extract the interesting samples #########################

interesting = c('MGH26', 'MGH28', 'MGH29', 'MGH30', 'MGH31')
mgh26Data = subsetSamples(data, 'MGH26')
mgh28Data = subsetSamples(data, 'MGH28')
mgh29Data = subsetSamples(data, 'MGH29')
mgh30Data = subsetSamples(data, 'MGH30')
mgh31Data = subsetSamples(data, 'MGH31')

############################################## GSVA #############################################

mgh26Signature = measureSignatures(mgh26Data, signatures)
mgh28Signature = measureSignatures(mgh28Data, signatures)
mgh29Signature = measureSignatures(mgh29Data, signatures)
mgh30Signature = measureSignatures(mgh30Data, signatures)
mgh31Signature = measureSignatures(mgh31Data, signatures)

############################################## Heatmaps #############################################
myPalette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

heatmap.2(t(mgh26Signature), cexRow=1.5, main="Enrichment of FACS marker signatures \n in Anoop et al", 
          keysize=1, trace="none", density.info="none", dendrogram="both", xlab="Samples", col=myPalette,
          offsetRow=c(1,1), margins=c(2,7.5), labRow=colnames(mgh26Signature), labCol=row.names(mgh26Signature))
