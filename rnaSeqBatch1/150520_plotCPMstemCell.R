library(ggplot2)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')

dat = read.delim("140213_normalisedLog_CPM_facs.txt")
genes = read.delim('140203_shortVSlong.txt')[,c(1,2)]

# Join on gene names
datJoin = merge.data.frame(genes, dat, by.x="ensembl_gene_id", by.y=0)

# List the genes which I want to plot expression scores for
gen = c("TUBB3", "PROM1", "GFAP", "NANOG", "NES", "OLIG2", "SOX2", "POU5F1", "LAMB1")

# Subset
datSub = datJoin[datJoin$external_gene_id %in% gen,]
datSub = datSub[c(1:5,7,8),]
row.names(datSub) = datSub$external_gene_id
datSub = datSub[,c(3:8)]