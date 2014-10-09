# I aquired data from the TCGA data freeze which has combined all the data together from the various plaforms
setwd("~/Documents/public-datasets/TCGA/methylation/")
list.files()

calls = read.delim("TCGA_GBM_dnameth_calls_20120112_ver3.txt")
scores = read.delim("TCGA_GBM_dnameth_scores_20120112_ver3.txt")
categories = read.delim("TCGA_GBM_dnameth_call_categories.txt")
data = read.delim("DNA.methylation.k6.txt")