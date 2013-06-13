#A script to sequentically read in all the files in a driectory, combine to a large dataframe and merge on gene to a single dataframe

setwd('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/')

#a function to define files
map.func = function(masterFileList) { #takes a dataframe listing the files and returns a list object of individual paitent gene scores
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

#define a function that collapses the dataframe by merging gene names. Uses merge_by function. Must navigate to the directory where the files are.
collapse.dataframe = function(bigDataFrame) {
  results = Reduce('merge_by', bigDataFrame)
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
  geneNames = bigDataFrame[[1]]$gene.symbol
  final.results = cbind(geneNames, trimmed.results)
  row.names(final.results) = final.results$geneNames
  return(final.results)
}