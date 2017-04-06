library(ggplot2)
library(reshape)
library(plyr)

setwd('~/Data/ATAC_seq_digitoninBulk/')
source('~/Code/Rscripts/Templates/multiplot.R')

dat = read.csv('~/Code/Rscripts/Postdoc/ATAC_seq/170404_bulkDigitonin/170406_atacDigitoninBulk.csv', header = T, row.names=1)

p1 <- ggplot(dat, aes(factor(chr), Percent_mtDNA, fill = factor(Lysis_agent), label = factor(Sample))) +
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Lysis buffer",
                      breaks=c("IGEPAL", "Digitonin"), labels=c("IGEPAL", "Digitonin")) +
    ggtitle("Percent mtDNA ATAC-seq reads") + geom_text(check_overlap = TRUE, size = 3) +
    xlab("groups") + ylab("% mtDNA reads") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=16))
p1

p2 <- ggplot(dat, aes(factor(chr), ESTIMATED_LIBRARY_SIZE, fill = factor(Lysis_agent), label = factor(Sample))) +
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Lysis buffer",
                      breaks=c("IGEPAL", "Digitonin"), labels=c("IGEPAL", "Digitonin")) +
    ggtitle("Percent mtDNA ATAC-seq reads") + geom_text(check_overlap = TRUE, size = 3) +
    xlab("groups") + ylab("% mtDNA reads") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=16))
p2

p3 <- ggplot(dat, aes(factor(chr), PERCENT_DUPLICATION, fill = factor(Lysis_agent), label = factor(Sample))) +
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Lysis buffer",
                      breaks=c("IGEPAL", "Digitonin"), labels=c("IGEPAL", "Digitonin")) +
    ggtitle("Percent mtDNA ATAC-seq reads") + geom_text(check_overlap = TRUE, size = 3) +
    xlab("groups") + ylab("% mtDNA reads") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=16))
p3