setwd('~/Documents/Cell_biology/proliferation/Resazurin/140614_summary')
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

calcDMSOcontrol = function(dataFrame) {
    vehicle = dataFrame[dataFrame$treatment %in% 'DMSO',]
    tmz = dataFrame[dataFrame$treatment %in% 'TMZ',]  
    tmz$rep1 = tmz$rep1 / vehicle$rep1
    tmz$rep2 = tmz$rep2 / vehicle$rep2
    tmz$rep3 = tmz$rep3 / vehicle$rep3
    tmz$mean = rowMeans(tmz[,c(4:6)], na.rm=T)
    tmz$sd = apply(tmz[,c(4:6)], 1, sd, na.rm=T)
    return (tmz)
}

calcProlifNormalised = function(dataFrame) {
    negative = dataFrame[dataFrame$cd133 %in% 'CD133_neg',]
    positive = dataFrame[dataFrame$cd133 %in% 'CD133_pos',]
    positive$rep1 = positive$rep1 / negative$rep1
    positive$rep2 = positive$rep2 / negative$rep2
    positive$rep3 = positive$rep3 / negative$rep3
    positive$mean = rowMeans(positive[,c(4:6)], na.rm=T)
    positive$sd = apply(positive[,c(4:6)], 1, sd, na.rm=T)
    return (positive)
}

extractPosNegReplicates = function(dataFrame) {
    neg = dataFrame[dataFrame$cd133 %in% 'neg',c(1,2,3,9)]
    pos = dataFrame[dataFrame$cd133 %in% 'pos',c(1,2,3,9)]
    negMean = mean(neg$dmsoCorrected)
    negSD = sd(neg$dmsoCorrected)
    posMean = mean(pos$dmsoCorrected)
    posSD = sd(pos$dmsoCorrected)
    negSummary = c(negMean, negSD)
    posSummary = c(posMean, posSD)
    result = rbind(negSummary, posSummary)
    result = as.data.frame(result)
    origin = c('negative', 'positive')
    result = cbind(origin, result)
    colnames(result) = c('origin', 'mean', 'sd')
    return (result)
}

normaliseCD133 <- function (dataFrame) {
    cd133Neg = dataFrame[dataFrame$cd133status %in% 'CD133_neg',]
    cd133Pos = dataFrame[dataFrame$cd133status %in% 'CD133_pos',]
    cd133NegAv = mean(cd133Neg$mean)
    cd133NegSd = sd(cd133Neg$mean) / sqrt(length(cd133Neg$mean))
    cd133PosAv = mean(cd133Pos$mean)
    cd133PosSd = sd(cd133Pos$mean) / sqrt(length(cd133Pos$mean))
    cd133 = as.data.frame(rbind(c(cd133NegAv, cd133NegSd), c(cd133PosAv, cd133PosSd)))
    cd133$cd133 = c('negative', 'positive')
    return (cd133)
}

################################ IO and subsetting ############################################

data = read.delim('140614_day7meanSD.txt')

# Subset the data for what stain I used
dataCD133 = data[data$cd133status %in% c("CD133_neg", "CD133_pos"),]
dataMatched = dataCD133[!dataCD133$clone %in% c("030a", "034a"),]
dataDoubleStain = data[!data$cd133status %in% c("CD133_neg", "CD133_pos"),]

################################ Plot growth data ############################################

growthPlot7 = ggplot(data=dataMatched[dataMatched$treatment %in% 'DMSO',], 
                     aes(x=clone, y=mean, fill=cd133status)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing proliferation at day 7 \nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

################################ Plot Temozolomide data ############################################

day7TMZ = calcDMSOcontrol(dataMatched)

tmzPlot7 = ggplot(data=day7TMZ, aes(x=clone, y=mean, fill=cd133status)) + 
    scale_fill_manual(values=c("blue", "yellow")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Cell number relative to DMSO control") +
    ggtitle("Comparing temozolomide sensitivty at day 7\nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################################## Normalise proliferation to CD133 negative #################################################

day7GrowthNorm = calcProlifNormalised(dataMatched[dataMatched$treatment %in% 'DMSO',])
day7GrowthNorm$ID = paste(day7GrowthNorm$clone, day7GrowthNorm$cd133status)

day7GrowthNormP = ggplot(data=day7GrowthNorm, aes(x=clone, y=mean, fill=populationBias)) + 
    #scale_fill_manual(values=c("darkorange", "royalblue","cyan","magenta")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Cell number relative to matched CD133 negative") +
    ggtitle("Proliferation of CD133 cells at day 7") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

##### T.test proliferation ####
dataCD133 = dataMatched[dataMatched$treatment %in% 'DMSO',]
t.test(mean ~ cd133status, data=dataCD133, paired=T)
dataCD133 = normaliseCD133(dataMatched)
dataCD133[2,1] = dataCD133[2,1] / dataCD133[1,1]
dataCD133[1,1] = 1
dataCD133[1,2] = 0.0136
dataCD133[2,2] = 0.03200569

day7GrowthCD133 = ggplot(data=dataCD133, aes(x=cd133, y=V1, fill=cd133)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Cell number relative to matched CD133 negative") +
    ggtitle("Proliferation of CD133 cells at day 7") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    annotate("text", label="p = 0.29", x=2, y=1.25, size = 4) # Annotate the plot with the p-value
day7GrowthCD133
############################################## Normalise TMZ to CD133 negative #################################################

day7TMZNorm = calcProlifNormalised(day7TMZ)

day7TMZNormP = ggplot(data=day7TMZNorm, aes(x=clone, y=mean, fill=populationBias)) + 
    scale_fill_manual(values=c("springgreen4", "magenta","cyan")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Cell number relative to matched CD133 negative") +
    ggtitle("Comparing temozolomide sensitivty at day 7\nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

multiplot(day7GrowthNormP, growthPlot7, tmzPlot7, day7TMZNormP, cols=2)

#### T.test TMZ ####
t.test(mean ~ cd133status, data=day7TMZ, paired=T)
dataTMZ = normaliseCD133(day7TMZ)

day7TMZCD133 = ggplot(data=dataTMZ, aes(x=cd133, y=V1, fill=cd133)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Cell number relative to CD133 negative") +
    ggtitle("Comparing temozolomide sensitivty at day 7\nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    annotate("text", label="p = 0.36", x=2, y=1.25, size =4) # Annotate the plot with the p-value
day7TMZCD133

############################################## Read in the ELDA results #################################################
# Clear memory
rm(list=ls())


normaliseCD133 <- function (dataFrame) {
  cd133Neg = dataFrame[dataFrame$cd133status %in% 'CD133_neg',]
  cd133Pos = dataFrame[dataFrame$cd133status %in% 'CD133_pos',]
  cd133NegAv = mean(cd133Neg$Estimate)
  cd133NegSd = sd(cd133Neg$Estimate) / sqrt(length(cd133Neg$Estimate))
  cd133PosAv = mean(cd133Pos$Estimate)
  cd133PosSd = sd(cd133Pos$Estimate) / sqrt(length(cd133Pos$Estimate))
  cd133 = as.data.frame(rbind(c(cd133NegAv, cd133NegSd), c(cd133PosAv, cd133PosSd)))
  cd133$cd133 = c('negative', 'positive')
  return (cd133)
}

source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
setwd('~/Documents/Cell_biology/microscopy/ELDA/140615_summary/')
data = read.delim('140615_ELDA_summary.txt')
data$group = paste(data$clone, data$cd133status)
# Remove recurrents
data = data[!data$clone %in% c('030a', '034a')]

eldaRawP = ggplot(data=data, aes(x=group, y=Estimate, fill=clone)) + 
    #scale_fill_manual(values=c("aquamarine1", "darkgoldenrod1")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(0.9)) +
    xlab("") + ylab("Sphere efficiency") +
    ggtitle("Raw ELDA data") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

percentData = data
percentData[,4:6] = as.data.frame(data[,4:6] ^ -1 *100)

############################################ Adjust for multiple testing ###############################################

percentData$adjust = p.adjust(percentData$Test, method='fdr')
# Add a column with stars describing if a test is significant
percentData$star <- " "
percentData$star[percentData$adjust < .05]  = "*"
percentData$star[percentData$adjust < .01]  <- "**"
percentData$star[percentData$adjust < .001] <- "***"

eldaSig = ggplot(data=percentData, aes(x=group, y=Estimate, fill=clone)) + 
    # Define my own colors
    #scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(position=position_dodge(), stat="identity", color='black') +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.2, position=position_dodge(.9)) +
    xlab("Clone") + ylab("Percent sphere formation") +
    # scale_fill_hue(name="CD133")+#, Legend label, use darker colors
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Sphere forming efficiency of GICs at day 7") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    geom_text(aes(label=star), colour="black", vjust=-5, size=10) +
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

eldaCD133 = normaliseCD133(percentData)
# Conduct t-test
t.test(Estimate ~ cd133status, percentData, paired=T)
eldaCD133$star <- " "
eldaCD133[2,'star']  = "*"

eldaSumm = ggplot(eldaCD133, aes(x=cd133, y=V1, fill=cd133)) + 
    # Define my own colors
    scale_fill_manual(values=c("firebrick1", "mediumseagreen")) +
    geom_bar(position=position_dodge(), stat="identity", color='black') +
    geom_errorbar(aes(ymin=V1-V2, ymax=V1+V2), width=.2, position=position_dodge(.9)) +
    xlab("CD133 status") + ylab("Percent sphere formation") +
    # scale_fill_hue(name="CD133")+#, Legend label, use darker colors
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Sphere forming efficiency of Glioma initiating cells") +
    scale_y_continuous(breaks=0:20*4) +
    # Setting vjust to a negative number moves the asterix up a little bit to make the graph prettier
    geom_text(aes(label=star), colour="black", vjust=-2.25, size=10) +
    annotate("text", label="p = 0.045", x=2, y=15, size =5) + # Annotate the plot with the p-value
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(eldaRawP, eldaSig, eldaSumm, eldaPercent, cols=2)