setwd('~/Documents/Cell_biology/metaFigure/140501_fig/')

# Just do day 7 for now

growth = read.delim('1400501_day7GrowthMatched.txt')
growth = growth[growth$treatment %in% 'growth',]
tmz = read.delim('140423_day7TMZprocessed.txt')
invasion = read.delim('140501_invasionSummary.txt')
invasion$clone = c('035_neg', '035_pos', '041_neg', '041_pos', '039_neg', '039_pos')
invasion$Patient = c('035', '035', '041', '041', '039', '039')
elda = read.delim('140501_ELDApercentData.txt')
elda = elda[!elda$Patient %in% c('030a', '034a'),]

write.table(growth, 'growth.txt', sep='\t')
write.table(tmz, 'tmz.txt', sep='\t')
write.table(elda, 'elda.txt', sep='\t')
write.table(invasion, 'invasion.txt', sep='\t')

# Mung the data in excel
###################################### Plot the heatmap ######################################
# Clear memory
rm(list=ls())
library(RColorBrewer)
setwd('~/Documents/Cell_biology/metaFigure/140501_fig/collationHeatMap/')

data = read.delim('140501_collatedCellBiol.txt')
heatData = as.matrix(data[,c(4,6,8,10)])

# Scale the matrix so the values look nice
heatData = scale(heatData)
heatDataT = t(heatData)

colnames(heatDataT) = c('020 -', '020 +', '035 -', '035 +', '039 -', '039 +','041 -', '041 +')
row.names(heatDataT) = c('Proliferation', 'Sphere formation', 'TMZ resistance', 'Invasion')

#Now call the heat map
corrdist = function(x) as.dist(1-cor(t(x)))
cc = brewer.pal(9, 'YlOrRd')

#pdf('140501_heatMap.pdf', paper='a4')
heatmap(heatDataT, col=cc, margins=c(8,5),cexRow=1.2, main='Summary of phenotypic properties of GIC panel', 
             xlab='Patient clones by CD133 status', ylab='', Rowv=NA)
#dev.off()

# I think clone 039 was swapped around FACS so rename it
heatDataT2 = heatDataT
colnames(heatDataT2) = c('020 -', '020 +', '035 -', '035 +', '039 +', '039 -','041 -', '041 +')

#pdf('140501_heatMapSwap039.pdf', paper='a4')
heatmap(heatDataT2, col=cc, margins=c(8,5),cexRow=1.2, main='If I swapped the #039 CD133 label', 
        xlab='Patient clones by CD133 status', ylab='', Rowv=NA)
#dev.off()

###################################### Now do the qPCR data ######################################
qPCR = read.delim('~/Documents/RNAdata/qPCRexpt/140331_cd133s/140505_matrix_dCT.txt', row.names=1)
colnames(qPCR) = c('011 -', '011 +', '030a -', '030a +', '035 -', '035 +', '041 -', '041 +', '020 -', '020 +', '030 -', '030 +')
qPCR = as.matrix(qPCR)

#pdf('140505_qPCR_heatMap.pdf', paper='a4')
heatmap(qPCR, col=cc, margins=c(8,5),cexRow=1.2, main='Summary of gene expression of GIC panel', 
        xlab='Patient clones by CD133 status', ylab='Gene', Rowv=NA, scale='row')
#dev.off()
