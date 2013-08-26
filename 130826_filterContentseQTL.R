#!/usr/bin/env Rscript
library(optparse)
#setwd('~/Documents/eQTL/130823_fullData/')

######################################################################################################################
#Take the overlapped snp file and gene expression matrix and filter each contents to make matrix eQTL more efficient
######################################################################################################################

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="An Rscript to that matches the column names in a mutation SNP like file with a gene expression matrix"),
  make_option(c("-m", "--mutationFile"), action="store", type='character', default='./somatic.filterMaf.mafToSNP.fixPatientNames.txt.matchPatientNames.txt',
              help="The mutation SNP-like file"),
  make_option(c("-g", "--geneExpressionMatrix"), action="store",type = 'character', default='./GBM.uncv2.mRNAseq_RSEM_normalized_log2.txt.matchPatientNames.txt',
              help="A tab delimited text file of the gene expression matrix"),
  make_option(c("-s", "--snpFiltered"), action="store", type='character', default='outputM.txt',
              help="The Filtered mutation file"),
  make_option(c("-c", "--cuttoff"), action="store",type = 'integer', default=10,
              help="The minimum percntage prevalence of the mutation in the cohort"),
  make_option(c("-o", "--outGene"), action="store", type='character', default='outputG.txt',
              help="The filtered gene expression matrix")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

countNAs = function(numericMatrix) {
  num = sum(is.na(numericMatrix))
  return (num)
}

mutations = read.delim(opt$mutationFile, row.names=1)
genes = read.delim(opt$geneExpressionMatrix, row.names=1)
genes$NAs = apply(genes, 1, countNAs)

#filter for the prevalence of a mutation in a patient
patientNumber = length(mutations)
mutFilter = subset(mutations, rowSums(mutations) >= patientNumber/opt$cuttoff)

#Filter the RNAseq data for expression. Used the RSEM log2 normalised data
geneFilter = subset(genes, genes$NAs <= patientNumber/4)
geneFilter = geneFilter[,c(1:patientNumber)]