getwd()
setwd('/Users/d.brown6/Documents/public-datasets/firehose/stddata__2013_12_10/GBM/20131210_dataReformatting/dataRearranging/')
list.files()

#################################### Data input and clustering visualisations ##############################################
data = read.delim('140110_agilentNoNulls.txt')
setwd('./correlations/')
dataM = as.matrix(data)

prom1 = dataM['PROM1',]
cd133Sig = cor(dataM, prom1, na.rm=T)