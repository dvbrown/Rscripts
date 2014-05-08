library(WGCNA)

correlateGeneWithGEM <- function (geneExpressionMatrix = dat, gene='PROM1') {
  # Calculate the correlation between PROM1 expression and all the genes in TCGA GBM
  # The geneExpressionMatrix should be a numeric matrix with genes as columns as rows and patients as rows
  # The gene should be a character string
  corrPval = corAndPvalue(x=geneExpressionMatrix[,gene], y=geneExpressionMatrix)
  
  # Extract the correlation and p-value from the returned list
  correlation = corrPval$cor
  # Measure the  coefficient of determination
  coeffDeter = correlation^2
  pVal = corrPval$p
  fdr = p.adjust(pVal, method='fdr')
  
  result = t(rbind(correlation, coeffDeter, pVal, fdr))
  colnames(result) = c('correlation', ' Rsquared', 'p-value', 'FDR')
  return (result)
}

plotCoexpression <- function (geneCorrelationMatrix, gene) {
  # The geneCorrelationMatrix should be the result of the correlateGeneWithGEM function
  # The gene should be a character string
  par(mfrow=c(2,2))
  hist(geneCorrelationMatrix[,1], main=paste(gene, 'coexpression'), breaks='FD', xlab='Weighted correlation values')
  hist(geneCorrelationMatrix[,2], main=paste(gene, 'coefficient of determination'), breaks='FD', xlab='Weighted correlation values')
  hist(geneCorrelationMatrix[,4], main=paste(gene,'Significance'), breaks='FD', xlab='FDR corrected p-values')
  plot((geneCorrelationMatrix[,2]), -log10(geneCorrelationMatrix[,4]), main=paste(gene, 'Correlation vs significance'), 
       ylab='-log10 FDR p-value', xlab='Absolute correlation', pch=20)
  par(mfrow=c(1,1))
}

makeSquareCoexpressionMatrix <- function (geneCorrelationMatrix, geneExpressionMatrix) {
  # Returns a square adjacency matrix containing the module of genes highly coexpression with the gene of interest
    geneNames = row.names(geneCorrelationMatrix)
  # Build the network adjacency
  # Use the top correlated genes with PROM1 and measure their correlation with the transcriptome
  longMatrix = adjacency(geneExpressionMatrix,  selectCols = geneNames, #for correlation networks only (see below); can be used to select genes whose adjacencies will be calculated. Should be either a numeric vector giving the indices of the genes to be used, or a boolean vector indicating which genes are to be used.
                             type = "unsigned", power = 6, corFnc = "cor", #corOptions = "use = 'p'",
                             distFnc = "dist", distOptions = "method = 'euclidean'")
  
  # Make the adjacency matrix square
  squareMatrix = (longMatrix[colnames(longMatrix),])
  return(squareMatrix)
}
