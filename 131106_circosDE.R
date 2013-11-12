setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/output/')

data = read.delim('131021_shortVSlong.txt')
bed = read.delim('~/Documents/public-datasets/annotationFiles/ensGeneID.bed')

data$genesID = row.names(data)

result = merge.data.frame(data, bed, by.x='ensembl_gene_id', by.y='Ensembl.Gene.ID')

output = result[,c(10, 11, 12, 4)]

write.table(output, '../../../../circos/131106_heatmapPilot.txt', row.names=F, quote=F)
