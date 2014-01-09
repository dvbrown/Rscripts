# This script performs the differential expression testing
load('140109_intoDataFramesWorkspace.RData')

library(limma)
setwd('~/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/')
list.files()

targets = readTargets('140109_targets.txt', sep='\t')