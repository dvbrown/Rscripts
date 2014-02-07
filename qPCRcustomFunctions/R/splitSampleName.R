splitSampleName <-
function(plateMap) {
  # The column with sample is the vector containing the sample names you wish to split
  splitted = colsplit.factor(plateMap[['sample']], split = " ", names = c('origin', 'gene'))
  result = cbind(plateMap, splitted)
  return (result)
}
