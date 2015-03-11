# Thge script wll analyse the sorted PDGCs that were sorted and then reanalysed
library(ggplot2)
source("~/Documents/Rscripts/multiplot.R")

subtractBaseline = function(dataFrame, baseline, value) {
    dataSorted = dataFrame[order(dataFrame[,'PDGC'], dataFrame[,'Subpopulation']),]
    baseline = dataSorted[dataSorted$Subpopulation %in% baseline,]
    treatment = dataSorted[dataSorted$Subpopulation %in% value,]
    subtract =  treatment[,c(3:6)] - baseline[,c(3:6)]
    result = cbind(treatment[,c(1,2)], subtract)
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

setwd("~/Documents/Rscripts/150310_sortReanalysis/")
dat = read.csv("rawData.csv", row.names=1)
freq = dat[,c(1:6)]
count = dat[,c(1,2,7:10)]

freqBase = rbind(subtractBaseline(freq, "mixed", "doubleNeg"), subtractBaseline(freq, "mixed", "CD44"),
                 subtractBaseline(freq, "mixed", "CD133"), subtractBaseline(freq, "mixed", "doublePos"))
freqBase = freqBase[order(freqBase[,'PDGC'], freqBase[,'Subpopulation']),]

countBase = rbind(subtractBaseline(count, "mixed", "doubleNeg"), subtractBaseline(count, "mixed", "CD44"),
                 subtractBaseline(count, "mixed", "CD133"), subtractBaseline(count, "mixed", "doublePos"))
countBase = countBase[order(countBase[,'PDGC'], countBase[,'Subpopulation']),]

# Convert from wide to long
freqLong = melt(freqBase, id.vars=c("PDGC", "Subpopulation"))
hist(freqLong$value, breaks="FD")

p035 = plotPDGC(freqLong, "MU035")
p020 = plotPDGC(freqLong, "MU020")
p039 = plotPDGC(freqLong, "MU039")

# Summarise based on biological replicates
bioRep = ddply(freqLong, c("Subpopulation", "variable"), summarise, meanChange = mean(value, na.rm=T), 
               sdChange = sd(value, na.rm=T),reps=length(value))
bioRep$seChange = bioRep$sdChange / (sqrt(bioRep$reps))

# Plot the biological replicates
bioP = ggplot(data=bioRep, aes(x=Subpopulation, y=meanChange, fill=variable)) +
            geom_bar(stat="identity", position=position_dodge(), colour="black") + 
            ggtitle("Reanalysis biological Replicates n = 3") +  
            scale_fill_manual(values=c("orange", "darkred", "darkblue", "forestgreen")) + 
            geom_errorbar(aes(ymin=meanChange-seChange, ymax=meanChange+seChange), width=.2, position=position_dodge(0.9)) +
            xlab("Subpopulation") + ylab("Percent difference \nrelative to mixed population") +
            theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(p020, p035, p039, bioP, cols=2)