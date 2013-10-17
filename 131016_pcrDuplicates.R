library(edgeR)
setwd('~/Documents/RNAdata/danBatch1/pcrDuplicates/')
source("~/Documents/Rscripts/120704-sortDataFrame.R")

dups = read.delim('countsGene_gtf_getDups.txt', skip=1)
notDups = read.delim('countsGene_gtf.txt', skip=1)

data = data.frame(dups$GIC_039_CAGATC.trim.bowtie.merged.sortS.sortP.getDup.bam, notDups$GIC_039_CAGATC.trim.bowtie.merged.sortS.sortP.rmDup.bam)
colnames(data) = c('DupsOnly', 'DupsRm')
row.names(data) = dups$Geneid

counts = DGEList(counts=data, group=c('dup', 'noDup'))
c = counts$counts
keep = rowSums(c >1) > 1
counts = counts[keep,]

c = as.data.frame(counts$counts)
cs = sort.dataframe(c, 1)