library(affyPLM)
library(limma)
library(sqldf)
library(annotate)
library(hgu133plus2.db)

source('~/Documents/Rscripts/120704-sortDataFrame.R')
setwd('~/Documents/public-datasets/rembrandt/rembrandt_GBM/')
allFiles = list.files()
allFiles 

# rawData = ReadAffy('00518392_U133P2.CEL', rm.mask=T, rm.outliers=T)
# exprData = exprs(rawData)
# rmaData = justRMA('00518392_U133P2.CEL')

#rmaData = just.rma(allFiles, rm.outliers=T, rm.mask=T, verbose=T, destructive=T)
#normalisedData = as.data.frame(exprs(rmaData))

# Open a connection to Test.sqlite database
# setwd('~/Documents/public-datasets/rembrandt/rembrandt_GBM/processedData/')
# db <- dbConnect(SQLite(), dbname="140624_rembrandtGBM.sqlite")
# dbWriteTable(conn = db, name = "rembrandtRMA", value = normalisedData, row.names = TRUE)
# dbDisconnect(db)

boxData = normalisedData[,150:200]
boxplot(boxData, col=rainbow(50), main='Rembrandt RMA normalised', par(las=2))

head(normalisedData)

################# Extract the probe ID mappings with the expression data #########
geneIDs = ls(hgu133plus2ENTREZID)
#get gene ID numbers from the annptation package allowing for multiple probes to match mulitple genes
geneSymbols <- as.character(unlist(lapply(mget(geneIDs,env=hgu133plus2SYMBOL),
                                          function (symbol) { return(paste(symbol,collapse="; ")) } )))
geneNames <- as.character(unlist(lapply(mget(geneIDs,env=hgu133plus2GENENAME),
                                        function (name) { return(paste(name,collapse="; ")) } )))
unigene <- as.character(unlist(lapply(mget(geneIDs,env=hgu133plus2UNIGENE),
                                      function (unigeneID) { return(paste(unigeneID,collapse="; ")) } )))
#strip the Hs from the start of unigene reference
unigene <- gsub("Hs\\.","",unigene)

#read the gene annotations into a dataframe for use in the topTable function of limma
genelist <- data.frame(GeneID=geneIDs,GeneSymbol=geneSymbols,GeneName=geneNames)
head(genelist)

################# Match the expression data with expression set object ########################## 
annotatedData = merge(normalisedData, genelist, by.x='row.names', by.y='GeneID')

################# Calculate row median then take the highest probe intensity #########
rowMed = apply(normalisedData, 1, median)
annotatedData$rowMed = rowMed
annotateSort = sort.dataframe(annotatedData, 232)

# Remove NA and duplicated geneIDs
annotateSort <- annotateSort[!is.na(annotateSort$GeneSymbol), ]
annotateSort <- annotateSort[!duplicated(annotateSort$GeneSymbol), ]
head(annotateSort)

dbWriteTable(conn = db, name = "rembrandtSummarised", value = annotateSort, row.names = TRUE)

write.table(annotateSort, './140624_rembrandtGBM_RMAsummarise.txt', sep='\t')

################# Write out the clinical Data ########################## 
clinicalData = read.csv('~/Documents/public-datasets/rembrandt/REMBRANDT_survivalinmonths.csv')
dbWriteTable(conn = db, name = "rembrandtClinical", value = clinicalData, row.names = TRUE)
dbDisconnect(db)