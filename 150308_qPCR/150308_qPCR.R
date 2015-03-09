library(reshape)
library(plyr)
library(ggplot2)
source('~/Documents/Rscripts/multiplot.R')
source('~/Documents/Rscripts/qPCRFunctions.R')

setwd("150308_qPCR/")

dat = read.delim("dat/150308_ctValuesMapped.txt", row.names=1)

# Summarise the data by replicate
datWide = ddply(dat, .(cDNA, Gene), summarise, rep1=Cp[1], rep2=Cp[2], meanRep = mean(Cp, na.rm=T), 
                sdRep = sd(Cp), reps=length(Cp))
# Remove the non interesting targets
datWide = datWide[c(1:112),]
View(datWide)

# Plot replicates
correlation = ggplot(data=datWide, aes(x=rep1, y=rep2, color=Gene)) + 
    geom_point(shape=19) + geom_smooth(method=lm, colour='red') +
    xlab("Replicate 1") + ylab("Replicate 2") + # Set axis labels
    ggtitle("Correlation of technical replicates") +  # Set title
    theme_bw(base_size=18)