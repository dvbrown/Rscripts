library(RColorBrewer)

#Using tools coverage the coverage of tags was calculated relative to the proximal promoter file from Yin et al 2013 Genome Biol
setwd('~/Documents/CREB/ChIPseqENCODE/insectBED1000/')

files = list.files(pattern='^wg')
 A549_1 = read.delim(files[1], header=F)
 A549_2= read.delim(files[2], header=F)
 Ecc1_1= read.delim(files[3], header=F)
 Ecc1_2= read.delim(files[4], header=F)
 Gm12878_1= read.delim(files[5], header=F)
 Gm12878_2= read.delim(files[6], header=F)
 hESC_1= read.delim(files[7], header=F)
 hESC_2= read.delim(files[8], header=F)
 Hepg2_1= read.delim(files[9], header=F)
 Hepg2_2= read.delim(files[10], header=F)
 K562_1= read.delim(files[11], header=F)
 K562_2= read.delim(files[12], header=F)

count = cbind(A549_1[,1:5], A549_2$V5, Ecc1_1$V5, Ecc1_2$V5, Gm12878_1$V5, Gm12878_2$V5, hESC_1$V5, hESC_2$V5,
              Hepg2_1$V5, Hepg2_2$V5, K562_1$V5, K562_2$V5)
colnames(count) = c('chromosome', 'start', 'end', 'name', 'A549_1','A549_2', 'Ecc1_1', 'Ecc1_2', 'Gm12878_1', 
                    'Gm12878_2', 'hESC_1', 'hESC_2', 'Hepg2_1', 'Hepg2_2', 'K562_1', 'K562_2')
count$name = paste(count$chromosome, count$start, sep=':')
count$sum = rowSums(count[,5:16])
count$mean = rowMeans(count[,5:16])
#subset the count table to remove intervals with no tags
editCount = subset(count, count$sum > 0)
sumCount = editCount[,c(1,2,3,4,17)]

#there are many duplicates for some reason
sumCount1 = unique(sumCount)
sumCount1$strand = '.'

heat = as.matrix(editCount[,5:16])
row.names(heat) = editCount$name
#add 1 so there are no 0s
heat1 = cor(heat)
plot(heat1)
cc = brewer.pal(12,'YlOrRd')
heatmap(heat1,cexRow=1, main='Correlation matrix CREB ChIP-seq', col=cc)

write.table(editCount, './130417_countIntersectPromoterInstance.txt', sep='\t', row.names=F)
write.table(sumCount1, './130417_interSectCountPromoter.txt', sep='\t', row.names=F)