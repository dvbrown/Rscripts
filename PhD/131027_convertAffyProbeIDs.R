library(hgu133plus2.db)
library(hgu133plus2cdf)
library(annotate)

# geneIDs = ls(hgu133plus2ENTREZID)
getGeneSymbol <- function (probeList) {
  # Return the geneSymbol fromhgu133plus2 probe ID
  geneIDs = as.character(probeList)
  #get gene ID numbers from the annptation package allowing for multiple probes to match mulitple genes
  geneSymbols <- as.character(unlist(lapply(mget(geneIDs,env=hgu133plus2SYMBOL),
                                            function (symbol) { return(paste(symbol,collapse="; ")) } )))
}


#read the gene annotations into a dataframe for use in the topTable function of limma
genelist <- data.frame(GeneID=geneIDs,GeneSymbol=geneSymbols,GeneName=geneNames)
rm(geneIDs, geneNames,geneSymbols,unigene)

files = list.files(pattern='*.csv')
f = lapply(files, read.csv, skip=1)

type1 = f[[1]][,1]
type1 = getGeneSymbol(type1)
write.csv(type1, 'Rchange/type1.csv')

type2 = f[[2]][,1]
type2 = getGeneSymbol(type2)
write.csv(type2, 'Rchange/type2.csv')

type1.1 = f[[3]][,1]
type1.1 = getGeneSymbol(type1.1)
write.csv(type1.1, 'Rchange/type1.1.csv')
type1.2 = f[[4]][,1]
type1.2 = getGeneSymbol(type1.2)
write.csv(type1.2, 'Rchange/type1.2.csv')
type1.3 = f[[5]][,1]
type1.3 = getGeneSymbol(type1.3)
write.csv(type1.3, 'Rchange/type1.3.csv')

type2.1 = f[[6]][,1]
type2.1 = getGeneSymbol(type2.1)
write.csv(type2.1, 'Rchange/type2.1.csv')
type2.2 = f[[7]][,1]
type2.2 = getGeneSymbol(type2.2)
write.csv(type2.2, 'Rchange/type2.2.csv')

