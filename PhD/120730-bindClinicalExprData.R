#Bind the expression matrix to the clinical data by rows on common columns
bindRowCommonCols = function(expressionMatrix, clinicalData) {
  #append the '.01' to the headers
  cols = colnames(clinicalData)
  cols = sub('$', '.01', cols)
  colnames(clinicalData) = cols
  #create a vector of the intersection of the dataframes
  union = intersect(expressionMatrix, clinicalData)
  #Use the intersection to subset the dataframes
  survival.cut = clinicalData[,c(union)]
  expression.cut = expressionMatrix[,c(union)]
  #append new rownames to the genelist
  genes = as.character(expressionMatrix$genes)
  genes.symbol = append(genes, 'barcode')
  genes.symbol = append(genes.symbol, 'survival')
  total = rbind(expression.cut, survival.cut)
  newTotal = cbind(genes.symbol, total)
  row.names(newTotal) = newTotal$genes1
  return (newTotal)
}

#the code needs debugging a list object is returned at some point.