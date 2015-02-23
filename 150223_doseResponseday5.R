library(plyr)
library(scales)
source('140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

summaryStats = function(dataFrame) {
    require(plyr)
    # Takes a dataframe with factor information and computes summary statistics based on levels of as least 2 factors
    # factor 1 and 2 are characters
    result <- ddply(dataFrame, c('Patient', 'ConcLog'), summarise,
                    N    = length(Value),
                    mean = mean(Value, na.rm=T),
                    sd   = sd(Value, na.rm=T),
                    se   = sd / sqrt(N) )
    return (result)
}

####################### IO day5 ###########################
setwd("~/Documents/Cell_biology/microenvironment/150216_doseResponse/")
list.files()
dasBos5 = read.delim("dasBos_day5.txt")
dasBos5$Group = paste(dasBos5$Patient, dasBos5$Drug, dasBos5$Conc, sep="_")
ruxIL5 = read.delim("ruxIL6.txt")
ruxIL5$Group = paste(ruxIL5$Patient, ruxIL5$Drug, ruxIL5$Conc, sep="_")

#### Summarise replicates ####
day5 = rbind(dasBos5, ruxIL5)
day5$ConcLog = log2(day5$Conc + 0.001)
day5$Value = day5$Value - 6000
rux = day5[day5$Drug %in% "Ruxolitinib",]
das = day5[day5$Drug %in% "Ruxolitinib",]
bos = day5[day5$Drug %in% "Bosutinib",]
il6 = day5[day5$Drug %in% "IL6",]

ruxM = summaryStats(rux)
dasM = summaryStats(das)
bosM = summaryStats(bos)
il6M = summaryStats(il6)

#### Plot the raw results ####

bosDay5 = ggplot(bosM, aes(x=ConcLog, y=mean, colour=Patient, group=Patient)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_point(stat="identity", position=position_dodge(), colour="black") + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Log 2 concentration (uM)") + ylab("Fluorescent intensity") +
    ggtitle("Dose response curve Bosutinib day 5") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
bosDay5

ruxDay5 = ggplot(ruxM, aes(x=ConcLog, y=mean, colour=Patient, group=Patient)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_point(stat="identity", position=position_dodge(), colour="black") + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Log2 Concentration (uM)") + ylab("Fluorescent intensity") +
    ggtitle("Dose response curve Ruxolitinib") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
ruxDay5

il6Day5 = ggplot(il6M, aes(x=ConcLog, y=mean, colour=Patient, group=Patient)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_point(stat="identity", position=position_dodge(), colour="black") + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Log2 Concentration (ng)") + ylab("Fluorescent intensity") +
    ggtitle("Dose response curve IL6") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
il6Day5

dasDay5 = ggplot(dasM, aes(x=ConcLog, y=mean, colour=Patient, group=Patient)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_point(stat="identity", position=position_dodge(), colour="black") + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Log2 Concentration (uM)") + ylab("Fluorescent intensity") +
    ggtitle("Dose response curve Dasatinib") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
dasDay5

multiplot(dasDay5, bosDay5, il6Day5, ruxDay5, cols=2)