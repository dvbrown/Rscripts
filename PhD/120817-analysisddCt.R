#analyze ddCt values
setwd('~/Documents/RNA data/qPCRexpt/120716-newPrimer/')
data = read.delim('120817-ddCt.txt')

#encapsulate the aggregation in a function
aggregateReplicates = function(dataFrame, dataColumn, sampleName) {
  data.mean = aggregate(dataColumn, by=list(sampleName), mean)
  data.sd = aggregate(dataColumn, by=list(sampleName), sd)
  final.data = merge(data.mean, data.sd, by.x='Group.1', by.y='Group.1')
  return (final.data)
}

#write a function to subset the data by the gene
selectGene = function(dataFrame, geneOfInterest) {
  gene = dataFrame[grep(geneOfInterest, row.names(dataFrame)),]
  return (gene)
}

b2m = aggregateReplicates(data, data$ddCt.B2M, data$sample) #aggregate the replicates to one value
gapdh = aggregateReplicates(data, data$ddCt.GAPDH, data$sample)
pkm2 = aggregateReplicates(data, data$ddCt.PKM, data$sample)

colnames(b2m) = c('sample', 'meanB2M', 'sdB2M')
colnames(gapdh) = c('sample', 'meanGAP', 'sdGAP')
colnames(pkm2) = c('sample', 'meanPKM', 'sdPKM')

poolData = merge(b2m, gapdh, by.x='sample', by.y='sample')
poolData = merge(poolData, pkm2, by.x='sample', by.y='sample')
row.names(poolData) = poolData$sample
rm(b2m, gapdh, pkm2)

b2m = selectGene(poolData, 'B2M')
creb = selectGene(poolData, 'CREB1')
fox = selectGene(poolData, 'FoxG1')
gapdh = selectGene(poolData, 'gapdh')
gfap = selectGene(poolData, 'gfap')
id1 = selectGene(poolData, 'ID1')
nanog = selectGene(poolData, 'Nanog')
notch = selectGene(poolData, 'Notch1')
oct4 = selectGene(poolData, 'Oct4')
pkm2 = selectGene(poolData, 'PKM2')
prom1 = selectGene(poolData, 'Prom1')
sox2 = selectGene(poolData, 'Sox2')
nes = selectGene(poolData, 'NES')
byGeneExp = as.list(b2m, creb, fox, gapdh, gfap, id1, nanog, notch, oct4, pkm2, prom1, sox2, nes)
geneList = c('(b2m', 'creb', 'fox', 'gapdh', 'gfap', 'id1', 'nanog', 'notch', 'oct4', 'pkm2', 'prom1', 'sox2', 'nes')