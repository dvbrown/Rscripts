library(ggplot2)
library(reshape)
library(plyr)

setwd('~/Documents/Presentations/2017/170302_GenteticsGenomics/')
source('~/Code/Rscripts/Templates/multiplot.R')

dat = read.csv('costSingle_cell.csv', row.names = 1)

plt = ggplot(data=dat, aes(x=Genes.per.cell, y=Cost.per.library, colour=Collection, 
                           label=factor(row.names(dat)))) +
    geom_point(shape=19, size=5) +
    ggtitle("") + geom_text(check_overlap = FALSE, size = 5, vjust=2) +
    scale_colour_manual(values=c("darkorange", "blue", "maroon", "forestgreen")) + 
    xlab("Genes detected per cell") + ylab("Cost per cell") +  
    scale_x_continuous(breaks=c(0,2000, 4000, 8000, 16000)) +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
plt