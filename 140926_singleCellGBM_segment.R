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
write.table(subsetGEM, './analysis/140926_gem4GSEA.txt', sep='\t')
write.table(subtype, './analysis/140926_phenotypes.txt', sep='\t', row.names=F)

##############################################  Vary the cutoffs for determining double positive and use GSEA ############################################## 
standard = cutOffDoubPos(sigData, 0.1, 0.1)
write.table(subtype, './analysis/GSEAcutoffs/phenotypes_standard.txt', sep='\t', row.names=F)

double0 = cutOffDoubPos(sigData, 0, 0)
write.table(double0, './analysis/GSEAcutoffs/phenotypes_double0.txt', sep='\t', row.names=F)

cd1330_cd4401 = cutOffDoubPos(sigData, 0, 0.1)
write.table(cd1330_cd4401, './analysis/GSEAcutoffs/phenotypes_0_01.txt', sep='\t', row.names=F)

cd133175_cd4425 = cutOffDoubPos(sigData, 0.175, 0.25)
length(cd133175_cd4425[cd133175_cd4425 == 'doublePositive'])
write.table(cd133175_cd4425, './analysis/GSEAcutoffs/phenotypes_0175_25.txt', sep='\t', row.names=F)

cd13302_cd4425 = cutOffDoubPos(sigData, 0.2, 0.25)
length(cd13302_cd4425[cd13302_cd4425 == 'doublePositive'])
write.table(cd13302_cd4425, './analysis/GSEAcutoffs/phenotypes_2_25.txt', sep='\t', row.names=F)