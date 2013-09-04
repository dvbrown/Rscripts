#!/usr/bin/Rscript
library(biomaRt)
library(org.Hs.eg.db)
library(optparse)

setwd('~/Documents/FredCSC/reformattedFiles/symbolConversion/')

geneNameToID <- function (inFile) {
  #Takes a matrix of gene Symbols as argument and returns the geneID mappings.
  #org.Hs.egALIAS2EG, org.Hs.egGENENAME, org.Hs.egSYMBOL2EG
  table = read.delim(inFile)
  IDvector = as.character(table[,1])
  holder = vector(mode = 'character', len=length(IDvector))
  i = 1
  #check that gene names match a gene ID
  for (gene in IDvector) {
    print(gene)
    x = try(as.character(mget(gene, org.Hs.egALIAS2EG)))
    if(class(x) == "try-error") {x = 'noMatch'} 
    else {holder[i] = x}
    i = i + 1
  }
  orfs = cbind(IDvector, holder)
  result = merge.data.frame(orfs, table, by.x='IDvector', by.y='From') 
  return (result)
}

option_list <- list(
  make_option(c("-e", "--explain"), action="store_true", default=FALSE,
              help="Takes the vector "),
  make_option(c("-i", "--inFile"), action="store",type = 'character', default='~/Documents/FredCSC/reformattedFiles/130829_gseaExpression.gct',
              help="A tab delimited text file of listing the geneID you wish to convert"),
  make_option(c("-o", "--outFile"), action="store", type='character', default='output.txt',
              help="The file you wish to output results as a tab delimited text file")
)
opt <- parse_args(OptionParser(option_list=option_list))

inFileZ = 'ziskin.txt'
inFileLgr = 'lgr5.txt'
inFileEph2 = 'ephB2.txt'
outFile = opt$outFile

out = geneNameToID(inFileZ)

result = merge.data.frame(y, data, by.x='IDvector', by.y='NAME')

write.table(result, outFile, sep='\t')