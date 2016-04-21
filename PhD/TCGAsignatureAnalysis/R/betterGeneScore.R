betterGeneScore <-
function(GeneExpressionMatrix, upRegGenes, downRegGenes) {
    #takes as input a z-score transformed gene expression matrix, and 2 dataFrames of upreg and downreg signature genes
    #subset the gene expression matrix with the signature
    upGenes = GeneExpressionMatrix[names(upRegGenes),]
    downGenes = GeneExpressionMatrix[names(downRegGenes),]
    upSig = colMeans(upGenes, na.rm=T)
    downSig = colMeans(downGenes, na.rm=T)
    score = upSig - downSig
    return (score)
}
