# This script will measure the proliferation over 3 days and also compare Resazurin and LDH
source('~/Documents/Rscripts/qPCRFunctions.R')

setwd('~/Documents/Cell_biology/proliferation/Resazurin/140318_testing')
data = read.delim('140324_data.txt', row.names=1)

# set the identifiers as factors
data$cells = as.factor(data$cells)
# Filter wells that had something in them
data = data[complete.cases(data$treatment),]

data$experiment = paste(data$cells, data$treatment)

data.split = split(data, data$experiment)

x = data.split[[1]]
x = rbind(x, c(x[1,1], x[1,2], rowMeans(x[,c(3:5)]), x[6,1], x[7,1], x[8,1]))
