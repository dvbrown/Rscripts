library(plyr)
source('140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

summaryStats = function(dataFrame, groups) {
    # The dataFrame with the readings
    # a vector of strings with the column names you wish to group by
    # Measurement is the numerical numbers you're interested in
    result = ddply(dataFrame, groups, summarise,
                   N    = sum(!is.na(fluorescence)),
                   mean = mean(fluorescence, na.rm=TRUE),
                   sd   = sd(fluorescence),
                   se   = sd / sqrt(N) )
    return (result)
}

############################################## IO ###############################################
setwd("~/Documents/Cell_biology/proliferation/Resazurin/141028_correct035/")
rawData = read.delim("pdgc035Reped.txt")
rawData = rawData[!is.na(rawData$subpop),]

dat = backgroundMeanSD(rawData, 7)
dat$cv = dat$sd / dat$mean * 100

dat = rawData
dat$mean = rowMeans(dat[c(4:6)])
dat$sd = apply(dat[c(4:6)], 1, sd)
dat$cv = dat$sd / dat$mean * 100

#### Plot the raw results ####
growthPlot7 = ggplot(dat[dat$treatment %in% 'DMSO',], 
                     aes(x=patient, y=mean, fill=subpop)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("PDAC") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
growthPlot7

temo = calcDMSOcontrol(dat)

tmzPlot7 = ggplot(data=temo, aes(x=patient, y=dmsoCorrected, fill=subpop)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("PDAC") + ylab("Cell number relative to \nDMSO control") +
    ggtitle("Temozolomide sensitivty at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))

multiplot(growthPlot7, tmzPlot7)