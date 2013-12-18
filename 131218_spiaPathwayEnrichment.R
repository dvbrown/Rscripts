# Measuring pathway enrichment using SPIA in RNA-seq batch1
library(SPIA)
source('~/Documents/Rscripts/131218_ensemblToEnterezConversion.R')
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/')

# Import the full dataset
data = read.delim('131021_shortVSlong.txt', row.names=2)
data = data[,c(2:8)]
data$ensembl = row.names(data)

# Add the enterez IDs
ensemblEnterezMap = ensembl2enterezConvert(data)
rnaseq.data.ent <- cbind(entid = ensemblEnterezMap$ENTREZID, data[ensemblEnterezMap$ENSEMBL,])
