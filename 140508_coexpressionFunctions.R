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
  pVal = correlation$p
  fdr = p.adjust(pVal, method='fdr')
  
  result = t(rbind(correlation, coeffDeter, pVal, fdr))
  colnames(result) = c('correlation', ' Rsquared', 'p-value', 'FDR')
  return (result)
}