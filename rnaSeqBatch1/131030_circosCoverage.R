setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/output/')

data = read.delim('131021_normalisedCPM.txt')
bed = read.delim('~/Documents/public-datasets/annotationFiles/ensGeneID.bed')

data$genesID = row.names(data)

result = merge.data.frame(data, bed, by.x='genesID', by.y='Ensembl.Gene.ID')

output = cbind(as.character(result$Chromosome.Name), result$Gene.Start..bp., result$Gene.End..bp.)

longA = as.character(rowMeans(result[,c(2,3,4)]))
shortA = as.character(rowMeans(result[,c(5,6,7)]))

outputL = cbind(output, longA)
outputS = cbind(output, shortA)

write.table(outputL, '../../../../circos/131030_rawCoverageAvLong.txt', row.names=F, quote=F)
write.table(outputS, '../../../../circos/131030_rawCoverageAvShort.txt', row.names=F, quote=F)