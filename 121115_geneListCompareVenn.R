library(VennDiagram)
setwd('~/Documents/RNAdata/RNAseqAnalysis/printFiles/121119_labMeeting/')
data = read.delim('121115_compareStudies.txt', stringsAsFactors=F)
data = data[,c(1,3,4,5,6)]
data1 = apply(data, c(1,2), gsub, '', NA)
dataList = list(data[,1],data[,2],data[,3],data[,4],data[,5])

venn.diagram(dataList[c(1,4,5)[]], filename='venn.tif', main='Overlap between CD133- and CD133+ gene expression studies', na='remove',
                 category.names=c('Yan 2011','Lottaz 2010','clone #035'),cat.pos=c(-20, 0, 20))

venn.diagram(dataList[c(2,3,5)[]], filename='venn2.tif', main='Overlap between CD133- and CD133+ gene expression studies', na='remove',
             category.names=c('Garcia 2010','Beier 2007','clone #035'),cat.pos=c(-20, 0, 20))

venn.diagram(dataList[c(5,4,1)[]], filename='venn3.tif', main='Overlap between CD133- and CD133+ gene expression studies', na='none',
             category.names=c('clone #035','Lottaz 2010','Yan 2011'),cat.pos=c(-20, 0, 20))