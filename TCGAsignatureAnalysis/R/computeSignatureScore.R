computeSignatureScore <-
function(zTransformedGem, geneSignature) {
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
