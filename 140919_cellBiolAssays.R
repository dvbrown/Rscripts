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
growthData = summaryStats(rawData, c('patient', 'subpop', 'treatment'))

############################################## IO ###############################################
setwd("~/Documents/Cell_biology/proliferation/Resazurin/")
pat004Data = read.delim("140821_004/pdac004Day7Label.txt", row.names=1)
pat041Data = read.delim("140919_gpsc041/pdac041Day7Label.txt", row.names=1)

rawData = rbind(pat004Data, pat041Data)
# Subtract background and take mean and SD
growthData = backgroundMeanSD(rawData)


#### Plot the raw results ####
growthPlot7 = ggplot(rawData=rawData[rawData$treatment %in% 'DMSO',], 
                     aes(x=patient, y=mean, fill=subpop)) + 
    scale_fill_manual(values=c("darkorange", "royalblue", "forestgreen", "red")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing growth at day 7 by \nmarker status") +  # Set title
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
growthPlot7