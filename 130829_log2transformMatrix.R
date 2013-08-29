#!/usr/bin/env Rscript
library(optparse)
#setwd('~/Documents/eQTL/130823_fullData/')

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help=To take an RNA-seq file and return the log2 version
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="To take an RNA-seq file and return the log2 version"),
  make_option(c("-g", "--geneExpressionMatrix"), action="store",type = 'character', default='~/Documents/eQTL/130829_fullData/illuminahiseq_RSEM_genes_normalized__data.data.txt',
              help="A tab delimited text file of the gene expression matrix"),
  make_option(c("-o", "--outGene"), action="store", type='character', default='outputG.txt',
              help="The log2 transformed gene expression matrix")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

geneMatrix = read.delim(opt$geneExpressionMatrix, row.names=1)

gaII = read.delim('~/Documents/eQTL/130829_fullData/GBM.uncv2.mRNAseq_RSEM_normalized_log2.txt')