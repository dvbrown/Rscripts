# A script to analyse the FACS experiment that I have done
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

setwd("~/Documents/facsData/flowJo/fortessa/150302_stupp/data/")
list.files()

cellCycle = read.delim("150305_cellCycleMunged.txt")
co = c("PDGC","Subpopulation", "Fraction","FreqOfParent","G1","S","G2")
cellCycle = cellCycle[,co]
col = c("darkgreen", "blue", "yellow", "darkred")

# Convert from wide to long
cellCycleLong = melt(cellCycle, id.vars=c("PDGC", "Subpopulation", "Fraction"))
frequency = cellCycleLong[cellCycleLong$variable %in% "FreqOfParent",]
cellCycleLong = cellCycleLong[!cellCycleLong$variable %in% "FreqOfParent",]

MU035 = ggplot(data=cellCycleLong[cellCycleLong$PDGC %in% "MU035",], 
              aes(x=variable, y=value, fill=Subpopulation)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    ggtitle("PDGC MU035") +  # Set title+
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU039 = ggplot(data=cellCycleLong[cellCycleLong$PDGC %in% "MU039",], 
               aes(x=variable, y=value, fill=Subpopulation)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    ggtitle("PDGC MU039") +  # Set title+
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU020 = ggplot(data=cellCycleLong[cellCycleLong$PDGC %in% "MU020",], 
               aes(x=variable, y=value, fill=Subpopulation)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    ggtitle("PDGC MU020") +  # Set title+
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU004 = ggplot(data=cellCycleLong[cellCycleLong$PDGC %in% "MU004",], 
               aes(x=variable, y=value, fill=Subpopulation)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    ggtitle("PDGC MU004") +  # Set title+
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

multiplot(MU020, MU035, MU039, MU004, cols=2)