setwd('~/Documents/public-datasets/RNA-seq/anoop2014_singleCellGBM/')

source('~/Documents/Rscripts/FACSmarkerTCGA/140508_coexpressionFunctions.R')

callMarkerSubtype = function(signatureScore, CD133cutoff, CD44cutoff) {
    # Takes a dataframe containing the signature scores and adds a new column that calls FACS marker subtype
    signatureScore$subtype = ""
    signatureScore$subtype = ifelse(signatureScore[,"CD133"] > signatureScore[,"CD44"], "CD133", "CD44")
    # signatureScore = sort.dataframe(signatureScore, "subtype")
    signatureScore$subtype = as.factor(signatureScore$subtype)
    return (signatureScore)
}

setwd('~/Documents/public-datasets/RNA-seq/anoop2014_singleCellGBM/')

sigData = read.delim('140926_signatureScoresAllPat.txt')
geneData = read.delim('GSE57872_GBM_data_matrix.txt', row.names=1)

##############################################  extract the double positive patients ############################################## 
doubPos = sigData[sigData$CD133 > 0.1 & sigData$CD44 > 0.1,]
double = row.names(doubPos)

subtyped = callMarkerSubtype(sigData, 0, 0)
subtyped$subtype = as.character(subtyped$subtype)
subtyped[double,]$subtype = "doublePositive"
subtyped$subtype = as.factor(subtyped$subtype)
subsetCases = row.names(subtyped)

##############################################  Subset the GEM and write phenotpye labels ############################################## 
subsetGEM = geneData[,subsetCases]
colnames(subsetGEM)
row.names(subtyped)

subtype = subtyped$subtype
write.table(subsetGEM, './analysis/140926_gem4GSEA.rnk', sep='\t')
write.table(subtype, './analysis/140926_phenotypes.cls', sep='\t')
