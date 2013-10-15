setwd('~/Documents/RNAdata/danBatch1/bowtieGem/')
files = list.files(pattern='*.txt')

f = lapply(files, read.delim, header=FALSE)
df1 = cbind(f[[1]],f[[2]],f[[3]],f[[4]],f[[5]],f[[6]])
df = df1[,c(2,4,6,8,10,12)]
row.names(df) = df1[,1]
colnames(df) = c('GIC_011', 'GIC_020', 'GIC_034', 'GIC_035', 'GIC_039', 'GIC_041')

noFeatures = tail(df)

dataSub = data[c(1:80921),]
row.names(dataSub) = data[,1]
dataSub$X15 = as.numeric(dataSub$X15)

s = sum(dataSub$X15)
tail(data)

data1 = dataSub[dataSub$X15 > 0,]