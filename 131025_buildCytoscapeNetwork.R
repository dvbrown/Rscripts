setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/')

geneData = read.delim('GLMedgeR/output/131021_shortVSlongLiberalDE.txt', row.names=1)
networkData = read.delim('networks/tabdelimited.H69__C2Mdbfx.txt')
geneData$description = gsub('Source:HGNC Symbol;Acc', '', geneData$description)

nw = networkData[,c(1,2,15)]

output = merge(nw, geneData, by.x='X.node1', by.y='external_gene_id')