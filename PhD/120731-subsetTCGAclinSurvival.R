creb.list = read.delim('120727-combineList3sources-annotated.txt')
creb.list = creb.list[,c(1,2)]
creb.list = as.character(creb.list$gene.symbol)
creb.list = unique(creb.list)

data.m = merge.data.frame(creb.list, expr.sub, by.x='gene.symbol', by.y='GeneSymbol')

clinData = read.delim('~/Documents/public-datasets/TCGA/clinicalData/120731-survivalDataStructered.txt')
tcga.matrix = read.delim('~/Documents/public-datasets/TCGA/expressionData/Expression-Genes/UNC__AgilentG4502A_07_1/120730-Agil1ExpressClinicalV1.txt')
row.names(tcga.matrix) = tcga.matrix$gene.symbol

#subset the tcga expression matrix with the gene list
 subUsingGeneList = function(expressionMatrix, geneList) {
  tcga.sub = expressionMatrix[geneList,] #subset the dataframe
  row.names(creb.sub) = creb.sub$genes
  #tcga.sub = cbind(geneList, tcga.sub)
  return (tcga.sub)
}
creb.sub = creb.sub[,2:123]

#transpose the data sp genes are now the columns
t.frame = t(creb.sub)
sum.stat = summary(t.frame)