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
das = day5[day5$Drug %in% "Dasatinib",]
bos = day5[day5$Drug %in% "Bosutinib",]
il6 = day5[day5$Drug %in% "IL6",]

ruxM = summaryStats(rux)
dasM = summaryStats(das)
bosM = summaryStats(bos)
il6M = summaryStats(il6)

#### Plot the raw results ####

bosDay5 = ggplot(data=bosM, aes(x = Conc, y = mean, group=Patient, colour = Patient)) + 
    geom_errorbar(width = 0.05, size = 0.75, aes(ymax = mean + sd, ymin = mean - sd, x = Conc)) +
    geom_point(alpha = 0.5, size = 5) +
    geom_line() +
    scale_x_log10("Log 2 concentration (uM)") + scale_y_continuous("Fluorescent intensity") +
    ggtitle("Dose response curve Bosatinib (Src, Abl) day 5") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                   text = element_text(size=22))

ruxDay5 = ggplot(data=ruxM, aes(x = Conc, y = mean, group=Patient, colour = Patient)) + 
    geom_errorbar(width = 0.05, size = 0.75, aes(ymax = mean + sd, ymin = mean - sd, x = Conc)) +
    geom_point(alpha = 0.5, size = 5) +
    geom_line() +
    scale_x_log10("Log 2 concentration (uM)") + scale_y_continuous("Fluorescent intensity") +
    ggtitle("Dose response curve Ruxitinib (Jak) day 5") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                   text = element_text(size=22))

il6Day5 = ggplot(data=il6M, aes(x = Conc, y = mean, group=Patient, colour = Patient)) + 
    geom_errorbar(width = 0.05, size = 0.75, aes(ymax = mean + sd, ymin = mean - sd, x = Conc)) +
    geom_point(alpha = 0.5, size = 5) +
    geom_line() +
    scale_x_log10("Log 2 concentration (pg/mL)") + scale_y_continuous("Fluorescent intensity") +
    ggtitle("Dose response curve IL-6 day 5") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                   text = element_text(size=22))

dasDay5 = ggplot(data=dasM, aes(x = Conc, y = mean, group=Patient, colour = Patient)) + 
    geom_errorbar(width = 0.05, size = 0.75, aes(ymax = mean + sd, ymin = mean - sd, x = Conc)) +
    geom_point(alpha = 0.5, size = 5) +
    geom_line() +
    scale_x_log10("Log 2 concentration (uM)") + scale_y_continuous("Fluorescent intensity") +
    ggtitle("Dose response curve Dasatinib (Src, Bcr-Abl, c-kit) day 5") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                   text = element_text(size=22))

multiplot(dasDay5, bosDay5, il6Day5, ruxDay5, cols=2)