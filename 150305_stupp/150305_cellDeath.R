# Analyse the dead cells
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

setwd("~/Documents/Rscripts/150305_stupp/dat/")
list.files()
col = c("darkgreen", "orange", "red", "yellow")
deadFreq = read.delim("150303_deadCells.txt")
deadCount = read.delim("150303_deadCellsCount.txt")
colnames(deadCount) = as.character(c("PDGC","Treatment","Dead_Count", "CD44-CD133-","CD44+","CD133+","CD44+CD133+"))
colnames(deadFreq) = as.character(c("PDGC","Treatment","Dead_Freq","CD44-CD133-","CD44+","CD133+","CD44+CD133+"))

deadCells = cbind(deadFreq[,c(1:3)], deadCount$Dead_Count)
colnames(deadCells) = c("PDGC", "Treatment", "Frequency", "Count")
deathLong = melt(deadCells, id.vars=c("PDGC", "Treatment"))

pCount = ggplot(data=deathLong[deathLong$variable %in% "Count",], 
                aes(x=PDGC, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Dead Cells (count)") + ggtitle("Absolute numbers") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

pFreq = ggplot(data=deathLong[deathLong$variable %in% "Frequency",], 
                aes(x=PDGC, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Dead Cells (percent of all cells)") + ggtitle("Percentage of Parent") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(pCount, pFreq, cols=1)

stuppTreated = deadFreq[deadFreq$Treatment %in% c("dmso","stupp"),c(1,2,4:7)]