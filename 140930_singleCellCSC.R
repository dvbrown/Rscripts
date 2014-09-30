# This is the same script as 140926_singleCellGBMSegment but for the stem cells

setwd('~/Documents/public-datasets/RNA-seq/anoop2014_singleCellGBM/analysis/')
source('~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')

callMarkerSubtype = function(signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}

cutOffDoubPos <- function (signatureScore, CD133cutoff, CD44cutoff) {
    doubPos = signatureScore[signatureScore$CD133 > CD133cutoff & signatureScore$CD44 > CD133cutoff,]
    double = row.names(doubPos)
    subtyped = callMarkerSubtype(signatureScore, 0, 0)
    
    subtyped$subtype = as.character(subtyped$subtype)
    subtyped[double,]$subtype = "doublePositive"
    subtyped$subtype = as.factor(subtyped$subtype)
    subsetCases = row.names(subtyped)
    subsetGEM = geneData[,subsetCases]
    colnames(subsetGEM)
    row.names(subtyped)
    result = subtyped$subtype
    return (result)
}

db <- dbConnect(SQLite(), dbname='~/Documents/public-datasets/RNA-seq/anoop2014_singleCellGBM/singleCellRNAseq.sqlite')
dbListTables(db)

sigData = dbReadTable(db, "inVitroSigScores", row.names=1)
geneData = dbReadTable(db, 'rawData', row.names=1)

standard = cutOffDoubPos(sigData, 0.1, 0.1)
write.table(standard, 'GSEAcutoffs/stemCellsInput/140930_standard.txt', sep='\t', row.names=F)

strict = cutOffDoubPos(sigData, 0.15, 0.15)
write.table(strict, 'GSEAcutoffs/stemCellsInput/140930_1515.txt', sep='\t', row.names=F)

moreStrict = cutOffDoubPos(sigData, 0.175, 0.175)
write.table(moreStrict, 'GSEAcutoffs/stemCellsInput/140930_175175.txt', sep='\t', row.names=F)

dbDisconnect(db) 