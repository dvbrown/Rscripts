library(ggplot2)
library(reshape)
library(plyr)

setwd('~/Data/ATAC-Seq_1/summary_info/')
source('~/Code/Rscripts/Templates/multiplot.R')
dat = read.csv('atac_seqsummary_metrics.csv')

dat = dat[c(6,9,14,11,13,10),]

dat$Percent_dups = dat$Percent_dups * 100
dat$percent_mtDNAreads = dat$percent_mtDNAreads * 100

# Remove proteinase RNase as it has a very high library size and a wierd looking bioanalyzer trace
#dat = dat[c(114,16,17),]

#### Barchart ####
color = rainbow(6)

lib = ggplot(data=dat, aes(x=Sample, y=EstimatedLibSize)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black", fill=color) + 
    ggtitle("Estimated library size") +
    scale_fill_manual(values=rainbow(17)) + 
    xlab("Library") + ylab("Estimated library size") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
lib

mt = ggplot(data=dat, aes(x=Sample, y=percent_mtDNAreads)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", fill=color) + 
  ggtitle("Percent mitochondiral reads") +
  scale_fill_manual(values=rainbow(17)) + 
  xlab("Library") + ylab("Percent mitochondiral reads") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
mt

dups = ggplot(data=dat, aes(x=Sample, y=Percent_dups)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", fill=color) + 
  ggtitle("Percent duplicated reads") +
  scale_fill_manual(values=rainbow(17)) + 
  xlab("Library") + ylab("Percent duplicated reads") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
dups

multiplot(lib, mt, dups, cols=2)