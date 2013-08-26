#!/usr/bin/env Rscript
library(optparse)
#setwd('~/Documents/eQTL/')

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="An Rscript to that matches the column names in a mutation SNP like file with a gene expression matrix")
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="An Rscript to that matches the column names in a mutation SNP like file with a gene expression matrix"),
  make_option(c("-m", "--mutationFile"), action="store", type='character', default='./Matrix_eQTL_R/130822_mutFileTestGoodNames.txt',
              help="The mutation SNP-like file"),
  make_option(c("-g", "--geneExpressionMatrix"), action="store",type = 'character', default='./Matrix_eQTL_R/130822_geneFileTest.txt',
              help="A tab delimited text file of the gene expression matrix"),
  make_option(c("-s", "--snpFiltered"), action="store", type='character', default='outputM.txt',
              help="The Filtered mutation file"),
  make_option(c("-o", "--outGene"), action="store", type='character', default='outputG.txt',
              help="The filtered gene expression matrix")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

mutations = read.delim(opt$mutationFile, row.names=1)
genes = read.delim(opt$geneExpressionMatrix, row.names=1)
intersectPatient = intersect(colnames(mutations), colnames(genes))

overlapGenes = genes[,intersectPatient]
overlapMuts = mutations[,intersectPatient]
row.names(overlapGenes) = row.names(genes)
row.names(overlapMuts) = row.names(mutations)

write.table(overlapGenes, opt$outGene, sep='\t', row.names=T)
write.table(overlapMuts, opt$snpFiltered, sep='\t',row.names=T)
