#A function to get the gene ontogeny of a list

setwd('~/Documents/CREB/Lembereger2008/')
data = read.delim('120914-overlapLembergerRamsay.txt', header=T)

#Takes a list of capitlised gene names as input
annotateListGoTerm = function(list) {
  setwd('~/Documents/public-datasets/annotationFiles/FLAT_FILES_072010/')
  ontogeny = read.delim('GENE_ONTOLOGY.txt')
  data.go = merge.data.frame(list, ontogeny, by.x='upper', by.y='symbol') #set the right column for the dataframe
  annotated = data.go[,c(1,5,6,7,9,11)]
  annotated$upper = as.character(annotated$upper)
  return (annotated)
}

#subset the tcga expression matrix with the gene list
subUsingGeneList = function(expressionMatrix, geneList) {
  tcga.sub = expressionMatrix[geneList,] #subset the dataframe
  row.names(tcga.sub) = tcga.sub$genes
  tcga.sub = unique(tcga.sub)
  return (tcga.sub)
}
subsetMatrix = subUsingGeneList(tcga.matrix, annotated$upper)
subsetMatrix.x = subsetMatrix.uni[c(1,2,4:45),2:92] #get rid of the NAs and the geneNames
write.table(subsetMatrix.x, '~/Documents/CREB/Lembereger2008/120914-tcgaExpressionSubset.txt', sep='\t')