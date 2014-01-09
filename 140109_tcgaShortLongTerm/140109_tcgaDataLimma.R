# This script performs the differential expression testing

library(limma)
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
list.files()

affy = read.delim('140109_affyMetrix.txt')
agilent = read.delim('140109_agilent.txt')
targets = readTargets('140109_targets.txt', sep='\t')