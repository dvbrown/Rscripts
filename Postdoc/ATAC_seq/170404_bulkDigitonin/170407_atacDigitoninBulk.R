library(ggplot2)
library(reshape)
library(plyr)

setwd('~/Data/ATAC_seq_digitoninBulk/')
source('~/Code/Rscripts/Templates/multiplot.R')

dat = read.csv('~/Code/Rscripts/Postdoc/ATAC_seq/170404_bulkDigitonin/170406_atacDigitoninBulk.csv', header = T, row.names=1)
dat$Percent_mtDNA = dat$Percent_mtDNA * 100

p1 <- ggplot(dat[!dat$Cell_line %in% "HCC38",], aes(factor(Lysis_agent), Percent_mtDNA, fill = factor(Cell_number))) +
    geom_dotplot(binaxis = "y", stackdir = "center") +#, position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Cell number",
                      breaks=c("5000", "50000"), labels=c("5k", "50k")) +
    ggtitle("") + xlab("Lysis buffer") + ylab("% mtDNA reads") +
    theme_bw(base_size=28) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=28))
p1

p2 <- ggplot(dat[!dat$Cell_line %in% "HCC38",], aes(factor(Lysis_agent), ESTIMATED_LIBRARY_SIZE, fill = factor(Cell_number))) +
    geom_dotplot(binaxis = "y", stackdir = "center") +#, position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Cell number",
                      breaks=c("5000", "50000"), labels=c("5k", "50k")) +
    ggtitle("") + xlab("Lysis buffer") + ylab("Estimated library size") +
    theme_bw(base_size=28) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=28))
p2

p3 <- ggplot(dat, aes(factor(chr), PERCENT_DUPLICATION, fill = factor(Lysis_agent), label = factor(Sample))) +
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Lysis buffer",
                      breaks=c("IGEPAL", "Digitonin"), labels=c("IGEPAL", "Digitonin")) +
    ggtitle("") + geom_text(check_overlap = TRUE, size = 3) +
    xlab("groups") + ylab("% PCR duplicates") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=16))
p3

#multiplot(p1,p2,p3, cols=1)