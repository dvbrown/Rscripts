library(ggplot2)
library(reshape)
library(plyr)

multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.delim(file=x,header=T, sep='\t', skip=6, nrows = 2)})
  return(datalist) }

setwd('~/Data/ATAC_seq_digitoninBulk/Insert_size/TextFile/')
file_names = list.files()
source('~/Code/Rscripts/Templates/multiplot.R')
percent = read.csv('~/Code/Rscripts/Postdoc/ATAC_seq/170404_bulkDigitonin/170404_atacDigitoninBulk.csv', row.names = 1)

dat = multmerge('~/Data/ATAC_seq_digitoninBulk/Insert_size/TextFile/Dedup/')
df <- do.call("rbind", dat)
df = df[c(1,3,5,7,9,11,13,15,17,19,21),]
row.names(df) = c("GM12878_Digitonin_0.03",'GM12878_5000_Kris_IGEPAL','GM12878_50,000_Dan_IGEPAL','GM12878_50,000_Sep_lysis','GM12878_Digitonin_0.01',
                  'GM12878_Digitonin_0.02','GM12878_IGEPAL_0.15','GM12878_IGEPAL_0.2','HCC38_IGEPAL_0.1','HCC38_Digitonin_0.02','GM12878_IGEPAL_0.1')

dat = cbind(percent, df)
write.csv(dat, '~/Code/Rscripts/Postdoc/ATAC_seq/170404_bulkDigitonin/170406_atacDigitoninBulk.csv')

p1 <- ggplot(dat, aes(factor(LIBRARY), PERCENT_DUPLICATION)) +
    geom_boxplot() + geom_jitter(aes(colour=Lysis_agent)) +
    ggtitle("Percent mtDNA ATAC-seq reads") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p1


p2 <- ggplot(dat, aes(factor(LIBRARY), ESTIMATED_LIBRARY_SIZE)) +
    geom_boxplot() + geom_jitter(aes(colour=Lysis_agent)) +
    ggtitle("Percent mtDNA ATAC-seq reads") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p2