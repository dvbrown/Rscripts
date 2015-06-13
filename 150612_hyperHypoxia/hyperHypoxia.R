library(ggplot2)
library(reshape)
library(plyr)
setwd("150612_hyperHypoxia/")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Treatment']),]
    baseline = dataSorted[dataSorted$Treatment %in% baseline,]
    treatment = dataSorted[dataSorted$Treatment %in% value,]
    subtract =  treatment[,c(3:6)] - baseline[,c(3:6)]
    result = cbind(treatment[,c(1,2)], subtract)
    return (result)
}

change = rbind(subtractBaseline(dat, "control", "h2o2"),
               subtractBaseline(dat, "control", "hypoxia"))


dat_all = read.csv("150612_hyperHypoxia.csv",row.names=1)
colnames(dat_all) = c("PDGC", "Treatment", "Viable", "CD44+/ CD133-",
                  "CD44+/ CD133+", "CD44-/ CD133+", "CD44-/ CD133-")

dat = dat_all[,c(1,2,4:7)]
dat = dat[!dat$Treatment %in% "isotype",]

# Convert from wide to long
mDat = melt(dat, id.vars = c("PDGC", "Treatment"))


# Biological replicates
bioRep = ddply(mDat, .(PDGC, Treatment, variable), summarise, meanValue = mean(value, na.rm=T), 
               sdDiff = sd(value, na.rm=T), reps=length(value))
bioRep$seDiff = bioRep$sdDiff / (bioRep$reps)

#### Barchart ####
# Plot the biological replicates
mu020 = ggplot(data=bioRep[bioRep$PDGC %in% "MU020",], 
              aes(x=variable, y=meanValue, fill=Treatment)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("Reanalysis biological Replicates n = 3") +  
    scale_fill_manual(values=c("orange", "darkred", "darkblue", "forestgreen")) + 
    geom_errorbar(aes(ymin=meanValue-seDiff, ymax=meanValue+seDiff), width=.2, position=position_dodge(0.9)) +
    xlab("Subpopulation") + ylab("Percent difference \nrelative to mixed population") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
mu020

### Boxplot ###
box = ggplot(data=mDat, aes(x=PDGC, y=value)) +
    geom_boxplot() + geom_point(aes(colour=Treatment), size =3, alpha=0.7,  position = position_jitter(w = 0.175)) +
    ggtitle("qPCR Summary") +
    xlab("Gene") + ylab("ddCt relative to CD133 negative") +
    scale_y_continuous(breaks = round(seq(-6, 6, by = 1),1)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(text = element_text(size=20))
box