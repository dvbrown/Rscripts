library(plyr)
library(scales)
source('140211_multiplotGgplot2.R')
source('~/Documents/Rscripts/cellBiologyAnalysisFunctions.R')

summaryStats = function(dataFrame) {
    require(plyr)
    # Takes a dataframe with factor information and computes summary statistics based on levels of as least 2 factors
    # factor 1 and 2 are characters
    result <- ddply(dataFrame, c('Patient', 'Conc'), summarise,
                    N    = length(Value),
                    mean = mean(Value, na.rm=T),
                    sd   = sd(Value, na.rm=T),
                    se   = sd / sqrt(N) )
    return (result)
}

####################### IO day3 ###########################
setwd("~/Documents/Cell_biology/microenvironment/150216_doseResponse/")
list.files()
dasBos3 = read.delim("DasBos_day3.txt")
dasBos3$Group = paste(dasBos3$Patient, dasBos3$Drug, dasBos3$Conc, sep="_")
# dasBos5 = read.delim("dasBos_day5.txt")
# dasBos7$Group = paste(dasBos7$Patient, dasBos7$Drug, dasBos7$Conc, sep="_")
ruxIL3 = read.delim("RuxIL_day3.txt")
ruxIL3$Group = paste(ruxIL3$Patient, ruxIL3$Drug, ruxIL3$Conc, sep="_")
# ruxIL5 = read.delim("ruxIL6.txt")
# ruxIL5$Group = paste(ruxIL5$Patient, ruxIL5$Drug, ruxIL5$Conc, sep="_")

#### Summarise replicates ####
day3 = rbind(dasBos3, ruxIL3)
#day3$Conc = log2(day3$Conc + 0.001)
rux = day3[day3$Drug %in% "Ruxolitinib",]
das = day3[day3$Drug %in% "Ruxolitinib",]
bos = day3[day3$Drug %in% "Bosutinib",]
il6 = day3[day3$Drug %in% "IL6",]

ruxM = summaryStats(rux)
dasM = summaryStats(das)
bosM = summaryStats(bos)
il6 = summaryStats(il6)

#### Plot the raw results ####
ruxDay3 = ggplot(ruxM, aes(x=Conc, y=mean, colour=Patient, group=Patient)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_point(stat="identity", position=position_dodge(), colour="black") + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Concentration (uM)") + ylab("Fluorescent intensity") +
    ggtitle("Dose response curve") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
ruxDay3

il6Day3 = ggplot(il6, aes(x=Conc, y=mean, colour=Patient, group=Patient)) + 
    scale_fill_manual(values=c("forestgreen", "royalblue", "darkorange", "red")) +
    # scale_fill_manual(values=c("black", "lightgrey", "darkgrey", "white")) +
    geom_point(stat="identity", position=position_dodge(), colour="black") + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Concentration (uM)") + ylab("Fluorescent intensity") +
    ggtitle("Dose response curve") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))
il6Day3