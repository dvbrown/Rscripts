setwd('~/Documents/public-datasets/TCGA/classficationSignature/')

da = read.delim('reFormatted_ClaNC840_centroids.txt')

d = da[,c(3:6)]
d$mode = apply(d, 1 ,median)

da$mode = d$mode

write.table(da, './131022_tryReformatThis.txt', sep='\t')