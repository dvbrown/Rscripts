setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/')

geneData = read.delim('GLMedgeR/output/131021_shortVSlongLiberalDE.txt', row.names=1)
networkData = read.delim('networks/tabdelimited.H69__C2Mdbfx.txt')
geneData$description = gsub('Source:HGNC Symbol;Acc', '', geneData$description)

nw = networkData[,c(1,2,15)]

output = merge(nw, geneData, by.x='X.node1', by.y='external_gene_id')
output = output[,c(1,2,3,4,5,6,9)]
colnames(output) = c('source', 'target', 'evScore', 'description', 'LFC', 'CPM', 'FDR')

write.table(output, './networks/131025_cytoscapeInput.txt', sep='\t', row.names=F, quote=F)

gd = geneData[,c(1,3,4,7)]
write.table(gd, './networks/131025_cytoscapeData.txt', sep='\t', row.names=F, quote=F)