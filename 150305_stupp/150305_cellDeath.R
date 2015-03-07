# Analyse the dead cells
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Treatment']),]
    baseline = dataSorted[dataSorted$Treatment %in% baseline,]
    treatment = dataSorted[dataSorted$Treatment %in% value,]
    subtract =  treatment[,c(3:7)] - baseline[,c(3:7)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

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

#### Calculate the pecentage difference ####
stupp = subtractBaseline(deadFreq, "dmso", "stupp")
tmz = subtractBaseline(deadFreq, "dmso", "tmz")
rad = subtractBaseline(deadFreq, "dmso", "rad")
dif = rbind(tmz, rad, stupp)[,c(1,2,4:7)]
diff = dif[order(dif$PDGC, dif$Treatment),]
diffLong = melt(diff, id.vars=c("PDGC", "Treatment"))

MU035 = ggplot(data=diffLong[diffLong$PDGC %in% "MU035",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU035") +  scale_y_continuous(limits=c(-20,20)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

MU020 = ggplot(data=diffLong[diffLong$PDGC %in% "MU020",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU020") +  scale_y_continuous(limits=c(-20,20)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

MU004 = ggplot(data=diffLong[diffLong$PDGC %in% "MU004",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU004") +  scale_y_continuous(limits=c(-20,20)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

MU039 = ggplot(data=diffLong[diffLong$PDGC %in% "MU039",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU039") +  scale_y_continuous(limits=c(-20,20)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(MU004, MU020, MU035, MU039, cols=2)