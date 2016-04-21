library(ggplot2)
library(plyr)
source('~/Documents/Rscripts/multiplot.R')

setwd('~/Documents/facsData/MUSE/Cell_Cycle/140819_cellCycle/')
data = read.csv('140819_cellCycle.csv')
data$GPSC = as.factor(data$GPSC)

events = ggplot(data=data, aes(x=Population, y=Events, fill=GPSC)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="GPSC number") +      # Set legend title
    #scale_y_continuous(breaks = round(seq(min(bindData$ddCt), max(bindData$ddCt), by = 0.5),0.5)) + # This modifies the scale of the y axis.
    xlab("Subpopulation") + ylab("Number of events") + # Set axis labels
    ggtitle("Number of events collected for Cell Cycle") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
events

g1 = as.data.frame(data$G1)
g1$phase = 'G1'
s = as.data.frame(data$S)
s$phase = 'S'
g2 = as.data.frame(data$G2)
g2$phase = 'G2'
colnames(g1) = c('percent', 'phase')
colnames(g2) = c('percent', 'phase')
colnames(s) = c('percent', 'phase')

gpscs = as.character(data$GPSC)
patients = as.character(c(gpscs, gpscs, gpscs))
cellcycle = rbind(g1, s, g2)
cellcycle = cbind(cellcycle, patients)

# write.table(cellcycle, '140824_mungThis.txt', sep='\t')

munged = read.delim('140824_munged.txt')
munged$gpsc = as.factor(munged$gpsc)

# Split this dataframe by patient
cycleByPatient <- split(munged, munged$gpsc)

# Stacked bar graph -- this is probably what you want
c035 = ggplot(data=cycleByPatient[[1]], aes(x=population, y=percent, fill=phase)) + geom_bar(stat="identity", colour="black") +
            scale_fill_manual(values=c("grey", "white", "black")) +
            xlab("Subpopulation") + ylab("Percent") + # Set axis labels
            ggtitle("GPSC patient #035 \ncell cycle distribution") +  # Set title+
            theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

c039 = ggplot(data=cycleByPatient[[2]], aes(x=population, y=percent, fill=phase)) + geom_bar(stat="identity", colour="black") +
    scale_fill_manual(values=c("grey", "white", "black")) +
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    ggtitle("GPSC patient #039 \ncell cycle distribution") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

c041 = ggplot(data=cycleByPatient[[3]], aes(x=population, y=percent, fill=phase)) + geom_bar(stat="identity", colour="black") +
    scale_fill_manual(values=c("grey", "white", "black")) +
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    ggtitle("GPSC patient #041 \ncell cycle distribution") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# multiplot(c035, c039, c041, events, cols=2)


cdata <- ddply(munged, c('phase', 'population'), summarise,
               N    = length(percent),
               mean = mean(percent),
               sd   = sd(percent),
               se   = sd / sqrt(N) )

repd = ggplot(data=cdata, aes(x=population, y=mean, fill=phase)) + geom_bar(stat="identity") +
            xlab("Subpopulation") + ylab("Percent") + # Set axis labels
            # geom_errorbar(position=position_dodge(width=0, height=3), width=.25, aes(ymin=mean-se, ymax=mean+se)) +
            ggtitle("Cell cycle distribution of 3 GPSC patients") +  # Set title+
            theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
repd

clusteredBar = ggplot(data=cdata, aes(x=phase, y=mean, fill=population)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    scale_fill_manual(values=c("white", "darkgrey", "grey", "black")) +
    xlab("Subpopulation") + ylab("Percent") + # Set axis labels
    geom_errorbar(position=position_dodge(0.9), width=.25, aes(ymin=mean-se, ymax=mean+se)) +
    ggtitle("Cell cycle distribution of 3 GPSC patients") +  # Set title+
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
clusteredBar

multiplot(c035, c039, c041, clusteredBar, cols=2)