library(sqldf)
library(GSVA)
library(gplots)
library(ggplot2)

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

plotHeatMap <- function (signatureScore, sampleName) {
    # Argument 1, a numeric matrix with the signature scores
    # Arg 2, a string representing the sampleName. Will be used for the file name
    fileName = paste(sampleName, 'heatMap.pdf', sep='_')
    pdf(fileName, width=11.69, height=8.27, useDingbats=FALSE)
    heatmap.2(t(signatureScore), cexRow=1.5, cexCol=0.6, main=sampleName, Rowv=colnames(signatureScore),
            keysize=1, trace="none", density.info="none", dendrogram="column", xlab="Samples", col=myPalette,
            offsetRow=c(1,1), margins=c(7,7), labRow=colnames(signatureScore), labCol=row.names(signatureScore))
  dev.off()
}

############################################## IO #############################################

dbListTables(db)
cd133Sig = dbReadTable(db, "cd133CuttOff")
cd44Sig = dbReadTable(db, "cd44CuttOff")
cd15 = dbReadTable(db, "cd15CuttOff")
signatures = list("CD133" = row.names(cd133Sig), "CD44" = row.names(cd44Sig)
                  ,'CD15' = row.names(cd15))
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
myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)

# plotHeatMap(mgh26Signature, 'MGH26')
# plotHeatMap(mgh28Signature, 'MGH28')
# plotHeatMap(mgh29Signature, 'MGH29')
# plotHeatMap(mgh30Signature, 'MGH30')
# plotHeatMap(mgh31Signature, 'MGH31')

############################################## Process the tumour bulk #############################################
bulk = data[,c('MGH26Tumor','MGH28Tumor', 'MGH29Tumor', 'MGH30Tumor', 'MGH31Tumor')]
bulkSignature = measureSignatures(bulk, signatures)
plotHeatMap(bulkSignature, 'Tumor bulk')

############################################## Plot all patients on an axis of CD133 and CD44 scores #############################################

# Stick the signature scores together
mgh26Signature = as.data.frame(mgh26Signature)
mgh28Signature = as.data.frame(mgh28Signature)
mgh29Signature = as.data.frame(mgh29Signature)
mgh30Signature = as.data.frame(mgh30Signature)
mgh31Signature = as.data.frame(mgh31Signature)
mgh26Signature$Patient = 'MGH26'
mgh28Signature$Patient = 'MGH28'
mgh29Signature$Patient = 'MGH29'
mgh30Signature$Patient = 'MGH30'
mgh31Signature$Patient = 'MGH31'

# Output
signatureScores = rbind(mgh26Signature, mgh28Signature, mgh29Signature, mgh30Signature, mgh31Signature)
# write.table(signatureScores, './140926_signatureScoresAllPat.txt', sep='\t')

cbPalette = c('cornflowerblue', 'darkgreen', 'darkred', 'magenta4', 'mediumblue')

# Draw scatterplot
ggplot(data=signatureScores, aes(x=CD133, y=CD44, color=Patient)) + 
    geom_point(shape=19, alpha=1) + geom_smooth(method=lm, colour='black') +
    scale_fill_manual(values=cbPalette) +
    xlab("CD133 signatures") + ylab("CD44 signatures") + # Set axis labels
    ggtitle("Anoop et al 2014 single cell RNAseq\nall patients by coexpression signature score") +  # Set title
    theme_bw(base_size=18)