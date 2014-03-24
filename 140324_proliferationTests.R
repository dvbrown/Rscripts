# This script will measure the proliferation over 3 days and also compare Resazurin and LDH

setwd('~/Documents/Cell_biology/proliferation/Resazurin/140318_testing')
data = read.delim('140324_data.txt', row.names=1)

# set the identifiers as factors
data$cells = as.factor(data$cells)