setwd('~/Documents/RNA data/qPCRexpt/')
#bind the well position to the plate annotation
annotation = read.csv('120716-newPrimer/120816-newPrimer.csv', header=FALSE)
ann = t(annotation)
wellPos = read.delim('120716-newPrimer/120817-384welllayouttranspose.txt', header=FALSE)
well = as.character(wellPos[1:198,])
map = cbind(well, ann)
row.names(map) = map[,1]

#read in the data
data = read.delim('./120716-newPrimer/120816-qPCR.txt', skip=1)
data = data[1:198,]
sub.data = data[,c(3,5,8)]
row.names(sub.data) = sub.data[,1]

#merge the 2 dataframes
full.data = merge(map, sub.data, by.x='row.names', by.y='row.names', sort=FALSE)

#clean up
rm(ann, annotation, well, wellPos)

write.table(full.data, './120716-newPrimer/120817-annotatedData.txt', sep='\t')