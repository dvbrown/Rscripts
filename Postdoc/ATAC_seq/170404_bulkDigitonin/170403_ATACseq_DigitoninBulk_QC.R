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

# Add the sample names
treatments = c("GM12878_Digitonin_0.03",'GM12878_5000_Kris_IGEPAL','GM12878_50,000_Dan_IGEPAL','GM12878_50,000_Sep_lysis','GM12878_Digitonin_0.01',
               'GM12878_Digitonin_0.02','GM12878_IGEPAL_0.15','GM12878_IGEPAL_0.2','HCC38_IGEPAL_0.1','HCC38_Digitonin_0.02','GM12878_IGEPAL_0.1')
percent = cbind(percent, treatments)
write.table(percent, '~/Code/Rscripts/Postdoc/170404_atacDigitoninBulk.csv', sep=",")
percent = read.csv('~/Code/Rscripts/Postdoc/ATAC_seq/170404_bulkDigitonin/170404_atacDigitoninBulk.csv', row.names = 1)

#### boxplot ####
color = rainbow(6)

p1 <- ggplot(percent, aes(factor(chr), Percent_mtDNA)) +
  geom_boxplot() + geom_jitter(aes(colour=Lysis_agent)) +
  ggtitle("Percent mtDNA ATAC-seq reads") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p1

p2 <- ggplot(percent, aes(factor(chr), Percent_mtDNA)) +
  geom_boxplot() + geom_jitter(aes(colour=Cell_line)) +
  ggtitle("Percent mtDNA ATAC-seq reads") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p2

multiplot(p1,p2, cols=2)