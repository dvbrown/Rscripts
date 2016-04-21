# Summary of qPCR data at this date
library(plyr)
setwd('~/Documents/RNAdata/qPCRexpt/140615_summary/')
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

summariseStatistics_ddCt <- function (dataFrame, groupVariableA='cd133', groupVariableB='gene.x') {
    # Enter arguments as character strings
    # Generate N, mean, sd and se statistics for a dataframe
    cData = ddply(dataFrame, c(groupVariableA, groupVariableB), summarise,
                  N    = sum(!is.na(ddCt)),
                  mean = mean(ddCt, na.rm=TRUE),
                  sd   = sd(ddCt, na.rm=TRUE),
                  se   = sd / sqrt(N) )    
    return (cData)
}
data = read.delim('140615_dCTddCTsummaryEdit.txt', row.names=1)
data$sample = paste(data$clone, data$cd133status, data$gene.x, sep='_')
data$origin = paste(data$clone, data$cd133status, sep='_')

deDupData = data[!duplicated(data$sample),]
dataPrimary = deDupData[!deDupData$clone %in% c('030a'),]
#write.table(dataPrimary, '140615_summarised_ddCT.txt', sep='\t')
dataCD133 = dataPrimary[dataPrimary$cd133status %in% "CD133_pos",]

interestingGenes = c('CREB1', 'LAMB1', 'MGMT', 'NANOG', 'NES', 'OCT4', 
                     'OLIG2', 'PROM1', 'SOX2', 'TUBB3', 'GFAP')

dataInteresting = dataCD133[dataCD133$gene.x %in% interestingGenes,]

ddCTplot = ggplot(data=dataInteresting, aes(x=gene.x, y=ddCt, fill=clone)) + 
    #scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Gene") + ylab("Expression normalised to CD133 negative") +
    ggtitle("Expression of stemness/ tumourigenicity markers in sorted GICs") +  # Set title
    coord_cartesian(ylim = c(0,12)) + geom_hline(yintercept=1, colour="red") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


################## Compute the average of the gene expressions ##############

dataSummary = summariseStatistics_ddCt(dataPrimary, groupVariableA='cd133status', groupVariableB='gene.x')
summaryInteresting = dataSummary[dataSummary$gene.x %in% interestingGenes,]

avPlot = ggplot(data=summaryInteresting[summaryInteresting$cd133status %in% "CD133_pos",], aes(x=gene.x, y=mean, fill=gene.x)) + 
    #scale_fill_manual(values=c("orangered1", "darkslateblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("Gene") + ylab("Expression normalised to CD133 negative") +
    ggtitle("Expression of stemness/ tumourigenicity markers in sorted GICs") +  # Set title
    coord_cartesian(ylim = c(0,20)) + geom_hline(yintercept=1, colour="red") +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


multiplot(ddCTplot, avPlot)
