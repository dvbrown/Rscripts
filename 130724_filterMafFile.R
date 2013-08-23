#!/usr/bin/env Rscript
library(optparse)
#setwd('~/Documents/eQTL/')

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="A script to filter the maf file into a vcf file? Some file that is good for the tools.")
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="An Rscript to that matches the column names in a mutation SNP like file with a gene expression matrix"),
  make_option(c("-m", "--mafFile"), action="store", type='character', default='./130823_fullData/step4_gbm_liftover.aggregated.capture.tcga.uuid.maf2.4.migrated.somatic.txt',
              help="The mutation TCGA file"),
  make_option(c("-o", "--output"), action="store", type='character', default='./130823_fullData/output.txt',
              help="The filtered maf file with fileds of interest")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

#########################################################################################
#A script to filter the maf file into a vcf file? Some file that is good for the tools.
#########################################################################################

mutations = read.delim(opt$mafFile, skip=4)

mutFilter = mutations[,c(1,9,16,26)]
#rm(mutations)

position = mutations[,c(1,5,6,7)]
position$Start_position = as.integer(position$Start_position)
position$End_position = as.integer(position$End_position)
position$Average_position = as.integer((position$Start_position + position$End_position)/2)
positionFinal = position[,c(1,2,5)]

mutFilter$Variant_Classification = as.character(mutFilter$Variant_Classification)
#remove those mutations classed as silent
mutClass = mutFilter[(mutFilter$Variant_Classification != 'Silent'),]
write.table(mutClass,opt$output,sep='\t',row.names=FALSE)