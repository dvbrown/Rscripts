#!/usr/bin/env Rscript

library(optparse)

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
  make_option(c("-e", "--explanation"), action="store_true", default=FALSE,
              help="A merging script"),
  make_option(c("-x", "--table1"), action="store", type='character', default=NULL,
              help="The first file to be merged in tab delimited format"),
  make_option(c("-y", "--table2"), action="store",type = 'character', default=NULL,
              help="The second table"),
  make_option(c("-o", "--outFile"), action="store", type='character', default='output.txt',
              help="The file you wish to output results as a tab delimited text file"),
  make_option(c("-a", "--colTable1"), action="store", type='integer', default=1,
              help="The column number of the first table to be merged by"),
  make_option(c("-b", "--colTable2"), action="store", type='integer', default=1,
              help="The column number of the second table to be merged by")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

##########################################################################
#A merging script
##########################################################################

