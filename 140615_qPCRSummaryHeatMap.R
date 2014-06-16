# Summary of qPCR data at this date
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
# Work on a t-test
summariseStatistics_ddCt(dataPrimary, groupVariableA='cd133status', groupVariableB='gene.x')

data = read.delim('140615_dCTddCTsummaryEdit.txt', row.names=1)
data$sample = paste(data$clone, data$cd133status, data$gene.x, sep='_')
data$origin = paste(data$clone, data$cd133status, sep='_')

deDupData = data[!duplicated(data$sample),]
dataPrimary = deDupData[!deDupData$clone %in% c('030a'),]
#write.table(dataPrimary, '140615_summarised_ddCT.txt', sep='\t')
dataCD133 = dataPrimary[dataPrimary$cd133status %in% "CD133_pos",]

ddCTplot = ggplot(data=dataCD133, aes(x=gene.x, y=ddCt, fill=clone)) + 
    #scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    xlab("Gene") + ylab("Expression normalised to CD133 negative") +
    ggtitle("Expression of stemness/ tumourigenicity markers in sorted GICs") +  # Set title
    #scale_y_continuous(breaks = round(seq(0, 20, by = 2),1)) + # This modifies the scale of the y axis.
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ddCTplot

################## Compute the average of the gene expressions ##############
dataAv = aggregate(ddCt ~ gene.x, data=dataCD133, mean)
dataSD = aggregate(ddCt ~ gene.x, data=dataCD133, sd)
dataAv = cbind(dataAv, dataSD$ddCt)

dataSummary = summariseStatistics_ddCt(dataPrimary, groupVariableA='cd133status', groupVariableB='gene.x')
colnames(dataAv) = c("gene", "mean", "sd")

avPlot = ggplot(data=dataSummary[dataSummary$cd133status %in% "CD133_pos",], aes(x=gene.x, y=mean, fill=gene.x)) + 
    #scale_fill_manual(values=c("orangered1", "darkslateblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(0.9)) +
    xlab("Gene") + ylab("Expression normalised to CD133 negative") +
    ggtitle("Expression of stemness/ tumourigenicity markers in sorted GICs") +  # Set title
    #scale_y_continuous(breaks = round(seq(0, 20, by = 2),1)) + # This modifies the scale of the y axis.
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
avPlot
