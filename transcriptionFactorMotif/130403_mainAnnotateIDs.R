#!/usr/bin/env Rscript
#Functions to convert between various ENSEMBL IDS
library(optparse)
option_list <- list(
    make_option(c("-e", "--explain"), action="store_true", default=FALSE,
                help="Takes the "),
    make_option(c("-k", "--key"), action="store",type='character', default='ensembl_transcript_id',
                help="The key you wish to search the ensembl database against. eg ensembl_transcript_id, external_gene_id(symbol) or entrezgene"),
    make_option(c("-i", "--inFile"), action="store",type = 'character', default='~/Documents/CREB/ChIPseqENCODE/nearestBED/test.bed.txt',
                help="A tab delimited text file of listing the geneSymbols of the promoters you wish to interrogate"),
    make_option(c("-o", "--outFile"), action="store", type='character', default='output.txt',
                help="The file you wish to output results as a tab delimited text file")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))
library(biomaRt)
library(org.Hs.eg.db)
source('/Users/d.brown6/Documents/Rscripts/transcriptionFactorMotif/130326_transcriptionFactorMatchingFunc.R')


bioMartKey = as.character(opt$key)
inFile = opt$inFile
outFile = opt$outFile
data = read.delim(inFile, header=F)

result = transcriptIDtoSymbol(data, bioMartKey)
write.table(result, outFile, row.names=F, sep='\t')