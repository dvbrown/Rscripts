#A script to sequentically read in all the files in a driectory, combine to a large dataframe and merge on gene to a single dataframe

setwd('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/')

#a function to define files
map.func = function(masterFileList) {
  dirs = masterFileList$File.Name
  dirs = as.character(dirs)
  entries = lapply(dirs, FUN=read.table ,sep='\t',quote='',header=T)
  entries
}

#define a custom merge function
merge_by = function(x, y) {
  pool.data = merge(x, y, by=c('gene.symbol'))
  pool.data
}

results = Reduce('merge_by', entries)

#generate a vector containing the paitent codes
expression.headers = aligentG4502A_07_1map$Sample

#select out of the dataframe only the numbers
trimmed.results = data[,grep("^value", colnames(data))]
#put in the correct column names
colnames(trimmed.results) = expression.headers
final.results = cbind('gene.name', expression.headers)