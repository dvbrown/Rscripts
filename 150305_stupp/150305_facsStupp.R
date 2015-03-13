# A script to analyse the FACS experiment that I have done
require(ggplot2)
library(reshape)
library(plyr)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Subpopulation']),]
    baseline = dataSorted[dataSorted$Subpopulation %in% baseline,]
    treatment = dataSorted[dataSorted$Subpopulation %in% value,]
    subtract =  treatment[,c(5:7)] - baseline[,c(5:7)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

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

bioRep = ddply(cellCycleLong, .(Subpopulation, variable), summarise, meanPhase = mean(value, na.rm=T), 
               sdPhase = sd(value, na.rm=T),reps=length(value))
bioRep$seDdCt = bioRep$sdPhase / (sqrt(bioRep$reps))

absolute = ggplot(data=bioRep, aes(x=Subpopulation, y=meanPhase, fill=variable)) +
            geom_bar(stat="identity", position=position_dodge(), colour="black") + 
            ggtitle("qPCR biological Replicates (n = 2 - 3)") +  scale_fill_manual(values=c("lightblue", "gold1", "orangered")) + 
            geom_errorbar(aes(ymin=meanPhase-seDdCt, ymax=meanPhase+seDdCt), width=.2, position=position_dodge(0.9)) +
            xlab("Gene") + ylab("ddCt") +
            theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Present the data relative to double negative ####

diff = rbind(subtractBaseline(cellCycle, "CD44-/CD133-", "CD44-/CD133-"),
             subtractBaseline(cellCycle, "CD44-/CD133-", "CD44+/CD133-"),
             subtractBaseline(cellCycle, "CD44-/CD133-", "CD44-/CD133+"),
             subtractBaseline(cellCycle, "CD44-/CD133-", "CD44+/CD133+"))

cellCycleLong = melt(diff, id.vars=c("PDGC", "Subpopulation"))

# Summarise the data by biological replicate
cellCycleLong = cellCycleLong[!cellCycleLong$PDGC %in% "MU039"]
bioRep = ddply(cellCycleLong, .(Subpopulation, variable), summarise, meanPhase = mean(value, na.rm=T), 
               sdPhase = sd(value, na.rm=T),reps=length(value))
bioRep$seDdCt = bioRep$sdPhase / (sqrt(bioRep$reps))

baseline = ggplot(data=bioRep, aes(x=variable, y=meanPhase, fill=Subpopulation)) +
            geom_bar(stat="identity", position=position_dodge(), colour="black") + 
            ggtitle("Cell cycle profile of PDGCs by subpopulation") +  scale_fill_manual(values=c("forestgreen", "lightblue", "gold1", "red")) + 
            geom_errorbar(aes(ymin=meanPhase-seDdCt, ymax=meanPhase+seDdCt), width=.2, position=position_dodge(0.9)) +
            xlab("Gene") + ylab("ddCt") +
            theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
baseline
