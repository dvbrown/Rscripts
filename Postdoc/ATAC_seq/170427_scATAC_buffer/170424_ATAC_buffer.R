library(ggplot2)
library(reshape)
library(plyr)

multmerge1 = function(mypath){
    filenames=list.files(path=mypath, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.delim(file=x,header=T, sep='\t', skip=6, nrows = 2)})
    return(datalist) }

multmerge2 = function(mypath){
    filenames=list.files(path=mypath, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.delim(file=x,header=F, sep='\t')})
    return(datalist) }

setwd('~/Data/ATAC/170420_ATACbuffer/LibMetrics/')
file_names = list.files()
source('~/Code/Rscripts/Templates/multiplot.R')

# Mung the library metrics
dat = multmerge1('~/Data/ATAC/170420_ATACbuffer/LibMetrics/')
df <- do.call("rbind", dat)
row.names(df) = list.files("./")

# Mung the mtDNA percentage
dat2 = multmerge2('~/Data/ATAC/170420_ATACbuffer/ChrsALignment/')
df2 <- do.call("cbind", dat2)
df3 = df2[,seq(2,73,2)]
row.names(df3) = df2[,1]
chrC = as.matrix(df3)
# Compute the sum
readSum = colSums(chrC)
percent = as.data.frame(chrC['chrM',] / readSum)
row.names(percent) = row.names(df)
colnames(percent) = "Percent_mtDNA"

# Bind the 2 together
datFrame = cbind(df, percent)
datFrame = datFrame[,c(2,4,5,7,10)]
rm(df, df2, df3, percent, chrC)
write.csv(datFrame, "~/Data/ATAC/170420_ATACbuffer/170426_atacMetrics.csv")

#### boxplot ####
color = rainbow(6)

p1 <- ggplot(datFrame, aes(factor(chr), Percent_mtDNA)) +
    geom_boxplot() + geom_jitter(aes(colour=Lysis_agent)) +
    ggtitle("Percent mtDNA ATAC-seq reads") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p1

p2 <- ggplot(datFrame, aes(factor(chr), Percent_mtDNA)) +
    geom_boxplot() + geom_jitter(aes(colour=Cell_line)) +
    ggtitle("Percent mtDNA ATAC-seq reads") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p2

multiplot(p1,p2, cols=2)