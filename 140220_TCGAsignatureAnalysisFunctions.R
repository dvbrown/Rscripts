library(survival)
library(limma)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

#objects for use in library
clinical = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/140109_clinicalDataTCGA.txt', row.names=1)
zSCore = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/140214_testSignaturesAgilent/140218_zTransformedTCGA.txt', sep='\t')

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
    status = survivalData$patient.followups.followup.vitalstatus
    status = gsub('alive', 0, status)
    status = gsub('dead', 1, status)
    status = as.numeric(status)
    deathFollow = (survivalData[,c(4,5)])
    censored = ifelse(deathFollow$days_to_death %in% NA, deathFollow$days_to_last_followup, deathFollow$days_to_death)
    return (censored)
}
betterGeneScore = function(GeneExpressionMatrix, upRegGenes, downRegGenes) {
    #takes as input a z-score transformed gene expression matrix, and 2 dataFrames of upreg and downreg signature genes
    #subset the gene expression matrix with the signature
    upGenes = GeneExpressionMatrix[names(upRegGenes),]
    downGenes = GeneExpressionMatrix[names(downRegGenes),]
    upSig = colMeans(upGenes, na.rm=T)
    downSig = colMeans(downGenes, na.rm=T)
    score = upSig - downSig
    return (score)
}

bindSignatureToSurvival = function(signatureScore, clinical) {
    #take the annotated clinical data and bind it a computed gene signature score for analysis
    clinical = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting//140109_clinicalDataTCGA.txt', row.names=1)
    
    clinical = clinical[,c(4,5)] #get just the survival and censorship. FIX DIS
    
    geneSet = as.data.frame(geneScore)
    set = intersect(row.names(geneSet), row.names(clinical)) #get the overlap
    data = cbind(geneSet[set,], clinical[set,]) #bind the 2 datasets together based on overlap
    colnames(data) = c('sigScore', 'survival', 'censorship')
    return(data)
}