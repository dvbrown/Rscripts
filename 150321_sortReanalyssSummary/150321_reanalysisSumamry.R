# Thge script wll summarise the 2 experiments I did that sorted and then reanalysed PDGCs
library(ggplot2)
library(plyr)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Subpopulation']),]
    baseline = dataSorted[dataSorted$Subpopulation %in% baseline,]
    treatment = dataSorted[dataSorted$Subpopulation %in% value,]
    subtract =  treatment[,c(4:7)] - baseline[,c(4:7)]
    result = cbind(treatment[,c(1:3)], subtract)
    return (result)
}

plotPDGC = function(dataFrame, pdgc="MU035") {
    # This function takes a dataFrane of some data and returns a ggplot object specific to the cell line supplied as a string
    extract = dataFrame[,"PDGC"] %in% pdgc
    datPlot = dataFrame[extract,]
    cols = c("orange", "darkred", "darkblue", "forestgreen")
    p = ggplot(data=datPlot, aes(x=Subpopulation, y=value, fill=variable)) +
        geom_bar(stat="identity", position=position_dodge(), colour="black") + 
        ggtitle(paste("Sorting an reanalysis after 7 days", pdgc, sep=" ")) +  scale_fill_manual(values=cols) + 
        xlab("Subpopulation") + ylab("Percent difference \nrelative to mixed population") +
        theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return (p) }

setwd("~/Documents/Rscripts/150321_sortReanalyssSummary/")
dat = read.csv("rawData.csv", row.names=1)
dat = dat[!dat$Subpopulation %in% "iso",]
expt1 = dat[dat$Date %in% "150310",]
expt2 = dat[dat$Date %in% "141106",]
# Remove the MU039 as it doesn't have all 4 subpops
expt2 = expt2[!expt2$PDGC %in% "MU039",]

#### Calculate and mung the data ####
freq1 = rbind(subtractBaseline(expt1, "mixed", "doubleNeg"), subtractBaseline(expt1, "mixed", "CD44"),
                 subtractBaseline(expt1, "mixed", "CD133"), subtractBaseline(expt1, "mixed", "doublePos"))

freq2 = rbind(subtractBaseline(expt2, "mixed", "doubleNeg"), subtractBaseline(expt2, "mixed", "CD44"),
              subtractBaseline(expt2, "mixed", "CD133"), subtractBaseline(expt2, "mixed", "doublePos"))
# Calculate MU039
mu039 = dat[dat$Date %in% "141106",]
baseline = mu039[2,c(4:7)]
mu039_dn = mu039["061114_039_DN.fcs", c(4:7)]
mu039_dn = mu039_dn - baseline
mu039_44 = mu039["061114_039_CD44.fcs", c(4:7)]
mu039_44 = mu039_44 - baseline
mu039Result = rbind(mu039_44, mu039_dn)
mu039Result = cbind(rbind(dat["061114_039_CD44.fcs", c(1:3)], dat["061114_039_DN.fcs", c(1:3)]),
                    mu039Result)

freq = rbind(freq1, freq2, mu039Result)
freq = freq[order(freq[,'PDGC'], freq[,'Subpopulation']),]

# convert to long
diff = melt(freq, id.vars=c("PDGC", "Subpopulation", "Date"))

######## Summarise by biological replicates ######## 
cols = c("gold", "red", "blue", "forestgreen")
bioRep = ddply(diff, .(PDGC, variable, Subpopulation), summarise, meanDiff = mean(value, na.rm=T), 
               sdDiff = sd(value, na.rm=T), reps=length(value))
bioRep$seDiff = bioRep$sdDiff / (bioRep$reps)
write.csv(bioRep, "150321_bioReps.csv")

ggplot(data=bioRep[bioRep$PDGC %in% "MU035",], 
       aes(x=Subpopulation, y=meanDiff, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("MU035 recurrent GBM n = 2") +  scale_fill_manual(values=cols) + 
    geom_errorbar(aes(ymin=meanDiff-seDiff, ymax=meanDiff+seDiff), width=.2, position=position_dodge(0.9)) +
    xlab("Sorted population") + ylab("% difference to mixed population") + scale_y_continuous(breaks = round(seq(-50, 50, by = 10),1)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Replicate across all PDGCs, EXCEPT for MU035 which is an outlier
pdgcRep = ddply(diff[!diff$PDGC %in% "MU035",], .(variable, Subpopulation), summarise, meanDiff = mean(value, na.rm=T), 
                sdDiff = sd(value, na.rm=T), reps=length(value))
pdgcRep$seDiff = pdgcRep$sdDiff / sqrt(pdgcRep$reps)

ggplot(data=pdgcRep, aes(x=Subpopulation, y=meanDiff, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    ggtitle("All PDGCs except MU035 \nn = 3 - 4") +  scale_fill_manual(values=cols) + 
    geom_errorbar(aes(ymin=meanDiff-seDiff, ymax=meanDiff+seDiff), width=.2, position=position_dodge(0.9)) +
    xlab("Sorted population") + ylab("% difference to mixed population") + scale_y_continuous(breaks = round(seq(-50, 50, by = 10),1)) +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
