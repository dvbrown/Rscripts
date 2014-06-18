getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
list.files()

#################################### Data input and clustering visualisations ##############################################
data = read.delim('140110_agilentNoNulls.txt')
setwd('./correlations/')
dataM = (as.matrix(data))

prom1 = dataM['PROM1',]
egfr = dataM['EGFR',]

dataM = t(data)

cd133 = cor(dataM)
save.image('1404089_correlationMatrix')