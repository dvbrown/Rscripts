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
  make_option(c("-c", "--cuttoffMut"), action="store",type = 'numeric', default=0.05,
              help="The minimum percentage prevalence of the mutation in the cohort. Default = 5%"),
  make_option(c("-n", "--cuttoffGene"), action="store",type = 'numeric', default=0.2,
              help="The maximum percentage of patients with NA expression of a gene in the cohort. Default = 20% "),
  make_option(c("-o", "--outGene"), action="store", type='character', default='outputG.txt',
              help="The filtered gene expression matrix"),
  make_option(c("-p", "--outMut"), action="store", type='character', default='outputM.txt',
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

#filter for the prevalence of a mutation in a patient. Remove genes with less n % mutation rate in GBM
patientNumber = as.numeric(length(mutations))
mutations$mutNum = rowSums(mutations)
mutFilter = subset(mutations, rowMeans(mutations) >= opt$cuttoffMut)

barplot(mutFilter$mutNum, names.arg=row.names(mutFilter), ylab='Number of patients with mutation',
        main='Frequency of mutations')

mutFilter = mutFilter[,c(1:patientNumber)]

#Filter the RNAseq data for expression. Remove genes with more than n NA values (ie not expressed). Used the RSEM log2 normalised data
geneFilter = subset(genes, genes$NAs <= patientNumber * opt$cuttoffGene)
geneFilter = geneFilter[,c(1:patientNumber)]

write.table(geneFilter, file=opt$outGene, sep='\t', row.names=T)
write.table(mutFilter, file=opt$outMut, sep='\t', row.names=T)