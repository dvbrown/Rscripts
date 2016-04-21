#A library to work with TCGA data
library(proxy)

#a function to define the files in a directory for reading in
map.fileList = function(masterFileList) {
  dirs = masterFileList$File.Name
  dirs = as.character(dirs)
  entries = lapply(dirs, FUN=read.table ,sep='\t',quote='',header=T)
  return(entries)
}

#define a custom merge function
merge_by = function(x, y) {
  pool.data = merge(x, y, by=c('gene.symbol'))
  return(pool.data)
}

#define a function that collapses the dataframe by merging gene names. Uses the above merge_by function
collapse.dataframe = function(bigDataFrame) {
  results = Reduce('merge_by', entries)
  #select out of the dataframe only the numbers
  trimmed.results = results[,grep("^value", colnames(results))]
  return(trimmed.results)
}

#A function to correctly label the columns with paitent barcodes
#Takes as input the map file containing sample names and the merged datafile of all paitents
labelSample = function(mapFile, trimmed.results, bigDataFrame) {
  #generate a vector containing the paitent codes
  expression.headers = as.character(mapFile$Sample)
  #put in the correct column names
  colnames(trimmed.results) = expression.headers
  #bind the gene names to the dataframe
  final.results = cbind(bigDataFrame$gene.symbol, trimmed.results)
  return(final.results)
}

#subset the tcga expression matrix with the supplied gene signature
subUsingGeneList = function(expressionMatrix, geneList) {
  tcga.sub = expressionMatrix[geneList,] #subset the dataframe
  row.names(tcga.sub) = tcga.sub$genes
  return (tcga.sub)
}

#Measure the similarity between a gene signature and the TCGA dataset
similarityMatrix = function(geneExpressionMatrix, GeneSignature) {
  index = 1
  distanceMatrix = geneExpressionMatrix[,,]
  while (index <= length(geneExpressionMatrix)) {
    distanceMatrix[,index] = simil(geneExpressionMatrix[,index], GeneSignature, pairwise=TRUE, by_rows=FALSE) #similarity measure
    index = index + 1
  }
  return (distanceMatrix)
}

#Wrap the the survival and similarity score data together for easy plotting
bindSurvivalSimilarity = function(distanceMatrix, survival) {
  final.data = cbind(distanceMatrix, survival)
  return (final.data)
}
