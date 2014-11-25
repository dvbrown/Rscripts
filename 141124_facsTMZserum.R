library(ggplot2)
library(plyr)

source("~/Documents/Rscripts/cellBiologyAnalysisFunctions.R")
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd("~/Documents/Cell_biology/proliferation/Resazurin/141124_tmzFcsFACS/")
color = c("chartreuse4", "skyblue2", "gold", "orangered1")
list.files()

growth = read.delim("141124_tmzFcsFACS.txt")
growth$sample = paste(growth[,"patient"], growth[,"treatment"], sep="_")

# Subtract background, take mean and sd
growth[,c(4:6)] = growth[,c(4:6)] - rowMeans(growth[17,c(4:6)])
growth$Mean = rowMeans(growth[,c(4:6)], na.rm=T)
growth$sd = apply(growth[,c(4:6)], 1, sd, na.rm=T)
growth$cv = growth$sd / growth$Mean * 100
growth = growth[c(1:16),]
write.table(growth, "141125_facsTMZserum.txt", sep='\t')

#Plot the raw results
growthPlot = ggplot(growth, aes(x=patient, y=Mean, fill=treatment)) + 
    scale_fill_manual(values=color) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("PDGC") + ylab("Fluorescent intensity") +
    ggtitle("Growth at day 7") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
growthPlot

growthSummary <- ddply(growth, "treatment", summarise,
                       N    = length(Mean), mean = mean(Mean),
                       sd   = sd(Mean), se   = sd / sqrt(N) )

growthSumPlot = ggplot(growthSummary, aes(x=treatment, y=mean, fill=treatment)) + 
    scale_fill_manual(values=color) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("PDGC") + ylab("Fluorescent intensity") +
    ggtitle("Average growth at day 7") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
growthSumPlot

multiplot(growthPlot, growthSumPlot, cols=1)