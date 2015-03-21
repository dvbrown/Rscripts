# A script to analyse the FACS experiment that I have done
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

setwd("~/Documents/facsData/flowJo/fortessa/150302_stupp/data/")
list.files()

dat = read.delim("150303_subpopTreatment.txt", row.names=1)
col = c("darkgreen", "orange", "red", "yellow")
colnames(dat) = as.character(c("PDGC","Treatment","CD44-CD133-","CD44+","CD133+","CD44+CD133+"))

# Convert from wide to long
datLong = melt(dat, id.vars=c("PDGC", "Treatment"))

MU035 = ggplot(data=datLong[datLong$PDGC %in% "MU035",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Percent of Parent") + ggtitle("PDGC MU035") +  # Set title+
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
subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Treatment']),]
    baseline = dataSorted[dataSorted$Treatment %in% baseline,]
    treatment = dataSorted[dataSorted$Treatment %in% value,]
    subtract =  treatment[,c(3:6)] - baseline[,c(3:6)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

stupp = subtractBaseline(dat, "dmso", "stupp")
tmz = subtractBaseline(dat, "dmso", "tmz")
rad = subtractBaseline(dat, "dmso", "rad")
dif = rbind(tmz, rad, stupp)
diff = dif[order(dif$PDGC, dif$Treatment),]
diffLong = melt(diff, id.vars=c("PDGC", "Treatment"))
col = c("darkorange", "darkred", "yellow")

MU035 = ggplot(data=diffLong[diffLong$PDGC %in% "MU035",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU035") +  #scale_y_continuous(limits=c(-50,50)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU020 = ggplot(data=diffLong[diffLong$PDGC %in% "MU020",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU020") +  #scale_y_continuous(limits=c(-50,50)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU039 = ggplot(data=diffLong[diffLong$PDGC %in% "MU039",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU039") +  #scale_y_continuous(limits=c(-50,50)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

MU004 = ggplot(data=diffLong[diffLong$PDGC %in% "MU004",], aes(x=variable, y=value, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=col) + 
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("MU004") +  #scale_y_continuous(limits=c(-50,50)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

multiplot(MU004, MU020, MU035, MU039, cols=2)

#### Summarise for some stats ####
require(plyr)
# Insert the second grouping variable here
diffSumry <- ddply(diffLong, c('variable', 'Treatment'), summarise,
                    N    = length(value), mean = mean(value),
                    sd   = sd(value), se   = sd / sqrt(N) )

colours = c("orange", "red", "yellow")
summryPlot = ggplot(data=diffSumry, aes(x=variable, y=mean, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=colours) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Subpopulation difference from vehicle") + 
    ggtitle("Plasticity of subtype with treatment \nn=4 PDGCs") +  #scale_y_continuous(limits=c(-50,50)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
summryPlot

# Remove MU004 and replot
diffSummry = ddply(diffLong[!diffLong$PDGC %in% "MU004",], c('variable', 'Treatment'), summarise,
                   N    = length(value), mean = mean(value),
                   sd   = sd(value), se   = sd / sqrt(N) )

sumryPlot = ggplot(data=diffSummry, aes(x=variable, y=mean, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +  
    scale_fill_manual(values=colours) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Difference to DMSO") + 
    ggtitle("Change of subtype n=3 PDGCs\n(no MU004)") +  #scale_y_continuous(limits=c(-50,50)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

multiplot(summryPlot, sumryPlot, cols=1)