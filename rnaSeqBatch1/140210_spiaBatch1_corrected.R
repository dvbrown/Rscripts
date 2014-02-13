# Measuring pathway enrichment using SPIA in RNA-seq batch1
library(SPIA)
library(pathview)
source('~/Documents/Rscripts/131218_ensemblToEnterezConversion.R')
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')

# Import the batch corrected cpm dataset
data = read.delim('140213_normalisedCPM_facs.txt')
#Import the log version
dataLog = read.delim('140213_normalisedLog_CPM_facs.txt')
#data = data[,c(2:8)]
#dataLog = dataLog[,c(2:8)]
data$ensembl = row.names(data)
dataLog$ensembl = row.names(dataLog)

# Add the enterez IDs
ensemblEnterezMap = ensembl2enterezConvert(data)
rnaseq.data.ent <- cbind(entid = ensemblEnterezMap$ENTREZID, data[ensemblEnterezMap$ENSEMBL,])

# Extract the differentially expressed genes
de.genes <- rnaseq.data.ent$logFC[data$FDR < 0.1]
names(de.genes) <- rnaseq.data.ent$entid[data$FDR < 0.1]
all.genes <- rnaseq.data.ent$entid

# Run spia
result.spia = spia(de=de.genes, all=all.genes, organism='hsa', nB=2000, plots=F)
result.spia$Name <- substr(spia.results$Name, 1, 25)
result.spia[1:20, -12]

# tA is the equivaent of fold hange pertubation in the pathway.
# pXXX is all the various FDR corrections

spiaNoNAs = result.spia[!is.na(result.spia$pPERT),]
plotP(spiaNoNAs, threshold=0.1)# x.lab='Enrichment score', ylab='Pertubation score', main='Disrupted pathways in short-term surviving GICs')

# output the results of the analysis
write.table(result.spia, './spia/140213_spiaResults.txt', sep='\t', row.names=F)
write.table(spiaNoNAs, './spia/140213_spiaNoNAsResults.txt', sep='\t', row.names=F)

# Visualise the KEGG pathway data in R
# Use the results of plotP to get the pathway ids
setwd('./spia/')



############### Repeat this analysis for log transformed batch corrected cpm and compare results #################
