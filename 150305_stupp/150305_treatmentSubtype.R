# A script to analyse the FACS experiment that I have done
require(ggplot2)
library(reshape)
source("~/Documents/Rscripts/multiplot.R")

setwd("~/Documents/facsData/flowJo/fortessa/150302_stupp/data/")
list.files()

dat = read.delim("150303_subpopTreatment.txt", row.names=1)
col = c("darkgreen", "orange", "red", "yellow")
colnames(dat) = as.character(c("PDGC","Treatment","CD44-CD133-","CD44+","CD133+","CD44+CD133+"))

# Convert from wide to long
datLong = melt(dat, id.vars=c("PDGC", "Treatment"))

MU035 = ggplot(data=datLong[datLong$PDGC %in% "MU035",], 
               aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent of Parent") + # Set axis labels
    ggtitle("PDGC MU035") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU020 = ggplot(data=datLong[datLong$PDGC %in% "MU020",], 
               aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent of Parent") + # Set axis labels
    ggtitle("PDGC MU020") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU039 = ggplot(data=datLong[datLong$PDGC %in% "MU039",], 
               aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent of Parent") + # Set axis labels
    ggtitle("PDGC MU039") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU004 = ggplot(data=datLong[datLong$PDGC %in% "MU004",], 
               aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=col) +
    xlab("Subpopulation") + ylab("Percent of Parent") + # Set axis labels
    ggtitle("PDGC MU004") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

multiplot(MU004, MU020, MU035, MU039, cols=2)

################# Calculate the percentage change #########################