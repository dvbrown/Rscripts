library(survival)
library(limma)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

extractDataClinical = function(clinicalData) { #A function to get useful clinical data and match rownames to expression
  sur = clinicalData[,c(2,7,8,10,13,19,29)]
  barcode = as.character(clinicalData$bcr_patient_barcode)
  barcode = gsub('$', '.01', barcode)
  sur = cbind(barcode, sur)
  row.names(sur) = sur$barcode
  return (sur)
}
#Arrange data so right censored paitents coded 0 with the days to death changed to days to last followup
censorData = function(survivalData) {
  status = survivalData$vital_status
  status = gsub('LIVING', 0, status)
  status = gsub('DECEASED', 1, status)
  status = as.numeric(status)
  deathFollow = (survivalData[,c(4,5)])
  censored = ifelse(deathFollow$days_to_death %in% NA, deathFollow$days_to_last_followup, deathFollow$days_to_death)
  return (censored)
}
betterGeneScore = function(GeneExpressionMatrix, upRegGenes, downRegGenes) {
  #takes as input a z-score transformed gene expression matrix, and 2 dataFrames of upreg and downreg signature genes
  upGenes = GeneExpressionMatrix[row.names(upRegGenes),]
  downGenes = GeneExpressionMatrix[row.names(downRegGenes),]
  upSig = colMeans(upGenes)
  downSig = colMeans(downGenes)
  score = upSig - downSig
  return (score)
}
#the tcga z Score data
tcga = read.delim('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/121109_AgilentPooledzTransform.txt', row.names=1)
#the gene signature
promoterCREB = read.delim('~/Documents/CREB/publicDataSets/Zhang2005-CREBMontminyDB/TableS1aPromPredictSortScore.txt')
promoterCREB = promoterCREB[1:200,c(1,2,5,16)] #the score column ranks the number of sites in the given promoter
koCREB = read.delim('~/Documents/CREB/publicDataSets/Lembereger2008/GSE8948_RAW/GPL8321/results/120816-Lemberger2008-normalVSCREBkoCutoff.txt')
koCREB = koCREB[,c(1,4,5,9,11)]
sortKO = sort.dataframe(koCREB, abs(koCREB$logFC), highFirst=T)
subsetTCGA = zScore[promoterCREB$Symbol,] #subset the tcga data with the gene list. This is the signature

upDown = (sort(apply(subsetTCGA, 1, median))) #get the median expression score to define up and down genes
#subset the gene score for upregulated and downregulated genes
up = upDown[upDown > 0] #upregulaed genes input for the geneSignature score function
down = upDown[upDown < 0]