library(ggplot2)
library(reshape)
library(plyr)

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.delim(file=x,header=F, sep='\t')})
  return(datalist) }

setwd('~/Data/ATAC_seq_digitoninBulk/Chromosome_aligned/')
file_names = list.files()
source('~/Code/Rscripts/Templates/multiplot.R')

dat = multmerge('~/Data/ATAC_seq_digitoninBulk/Chromosome_aligned/')
df <- do.call("cbind", dat)
names = rep(file_names, each = 2)

# Mung data into useful for for analysing reads on chromosomes
colnames(df) = names
row.names(df) = df[,1]
df = df[c(1:25),]
chromoCounts = df[,c(2,4,6,8,10,12,14,16,18,20,22)]
chrC = as.matrix(chromoCounts)

# Compute the sum
readSum = colSums(chrC)

percent = as.data.frame(chrC['MT',] / readSum)
colnames(percent) = ('Percent_mtDNA_reads')
percent$chr = 'mtDNA'

#### Barchart ####
color = rainbow(6)

p <- ggplot(percent, aes(factor(chr), Percent_mtDNA_reads)) +
  geom_boxplot() + geom_jitter() +
  ggtitle("Percent mtDNA ATAC-seq reads") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p

