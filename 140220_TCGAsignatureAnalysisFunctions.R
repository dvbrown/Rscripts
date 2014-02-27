library(survival)
library(limma)
source('~/Documents/Rscripts/120704-sortDataFrame.R')

#objects for use in library
clinical = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/140109_clinicalDataTCGA.txt', row.names=1)
zScore = read.delim('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/140214_testSignaturesAgilent/140218_zTransformedTCGA.txt', sep='\t', row.names=1)

#Arrange data so right censored paitents coded 0 with the days to death changed to days to last followup
censorData = function(survivalData) {
    # Takes patient dataframe
    status = survivalData[,'patient.followups.followup.vitalstatus']
    status = gsub('alive', 0, status)
    status = gsub('dead', 1, status)
    status = gsub('deceased', 1, status)
    # Check if the survival analysis can handle NAs for survival
    result = survivalData[,c(1,4,5)]
    censored = cbind(result, status)
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
computeSignatureScore = function(zTransformedGem, geneSignature) {
    # Takes a gene expression matrix that has been z transformed and a vector of gene names to return a signature score for each patient
    #Code that checks for membership of the gene signature in the zScore
    overlap = intersect(row.names(zTransformedGem), geneSignature)
    subsetTCGA = zTransformedGem[overlap,] #subset the tcga data with the gene list. This is the signature     
    upDown = (sort(apply(subsetTCGA, 1, median))) #get the median expression score to define up and down genes
    #subset the gene score for upregulated and downregulated genes
    up = upDown[upDown > 0] #upregulaed genes input for the geneSignature score function
    down = upDown[upDown < 0]
    #compute the geneScore
    geneScore = betterGeneScore(zTransformedGem, up, down)
    # Returns a vector with the scores
    return (geneScore)
}
bindSignatureToSurvival = function(signatureScore, clinical) {
    #take the annotated clinical data and bind it a computed gene signature score for analysis    
    clinical = clinical[,c(1,4)] #get just the survival and censorship. FIX DIS
    geneSet = as.data.frame(signatureScore)
    set = intersect(row.names(geneSet), row.names(clinical)) #get the overlap
    data = cbind(geneSet[set,], clinical[set,]) #bind the 2 datasets together based on overlap
    colnames(data) = c('sigScore', 'survival', 'censorship')
    return(data)
}
buildClassifier = function(signatureSurvivalFrame, percentileNum) {
    # Take a dataframe containing a signature score and censorship status and add the group membership eg 'high' or 'low'
    # Allows one to vary the percentile used as the classifier
    percent = quantile(signatureSurvivalFrame$sigScore, probs=percentileNum, names=T)
    signatureSurvivalFrame$percentile = ifelse(signatureSurvivalFrame$sigScore >= percent, 'high', 'low')
    return (signatureSurvivalFrame)
}

package.skeleton(name = 'TCGAsignatureAnalysis', path='~/Documents/Rscripts/', force=F) 