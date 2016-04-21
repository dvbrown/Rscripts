source('~/Documents/Rscripts/140211_multiplotGgplot2.R')
setwd('~/Documents/Cell_biology/microscopy/invasion/140606_4clones/')

invD7 = read.delim("140616_processedV1.txt")
invD7$clone = c('004','004','034','034','004','004','034','034','020','020',
               '028','028','020','020','028','028','020','020','028','028','020','020')
invD7$clone = as.factor(invD7$clone)

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data take the mean and sd
    dataFrame$mean = rowMeans(dataFrame[,c(5:7)], na.rm=T)
    dataFrame$sd = apply(dataFrame[,c(5:7)], 1, sd, na.rm=T)
    return (dataFrame)
}

normaliseMatrixCD133 = function(dataFrame) {
    # Normalise Matrix first
    noMatrix = dataFrame[dataFrame$matrix %in%  FALSE,]
    matrix = dataFrame[dataFrame$matrix %in% TRUE,]
    
    matrix$rep1 = matrix$rep1 / noMatrix$rep1
    matrix$rep2 = matrix$rep2 / noMatrix$rep2
    matrix$rep3 = matrix$rep3 / noMatrix$rep3
    matrix$mean = rowMeans(matrix[,c(5:7)], na.rm=T)
    matrix$sd = apply(matrix[,c(5:7)], 1, sd, na.rm=T)
    
    # Normalise CD133
    negative = matrix[matrix$cd133 %in% 'CD133_neg',]
    positive = matrix[matrix$cd133 %in% 'CD133_pos',]
    positive$rep1 = positive$rep1 / negative$rep1
    positive$rep2 = positive$rep2 / negative$rep2
    positive$rep3 = positive$rep3 / negative$rep3
    positive$mean = rowMeans(positive[,c(5:7)], na.rm=T)
    positive$sd = apply(positive[,c(5:7)], 1, sd, na.rm=T)
    # Return both dataframes
    result = list(matrix, positive)
    return (result)
}
############################################## Preprocessing #################################################


invD7 = invD7[!invD7$matrix %in% NA,]
invD7$sample = paste(invD7$clone, invD7$cd133, sep='_')

invD7Stats = backgroundMeanSD(invD7)
# Remove NAs
invD7Stats = invD7Stats[!invD7Stats$mean %in% NaN,]

invD7StatsP = ggplot(data=invD7Stats[!invD7Stats$clone %in% c('030a', '034a'),], aes(x=sample, y=mean, fill=matrix)) + 
    scale_fill_manual(values=c("skyblue3", "yellow")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere") +
    ggtitle("Invasive ability of GIC clones by CD133 status at day 7") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################################## Normalise for no matrix and CD133 #################################################

invD7Norm = normaliseMatrixCD133(invD7Stats)

# noMatrix = invD7Stats[invD7Stats$matrix %in%  FALSE,] This part is to calculate 020 neg for the heatmap
# noMatrix = noMatrix[c(1:7),]
# matrix = invD7Stats[invD7Stats$matrix %in% TRUE,]
# matrix$matNormalised = matrix$mean / noMatrix$mean
# matrix$matNormalisedSD = matrix$sd / noMatrix$sd

matchedCd133 = invD7Norm[[1]][invD7Norm[[1]]$cd133status %in% c('CD133_neg', "CD133_pos"),]
doubleStain = invD7Norm[[1]][!invD7Norm[[1]]$cd133status %in% c('CD133_neg', "CD133_pos"),]

# write.table(matchedCd133, "140616_matchedCD133Invasion.txt", sep='\t')
# write.table(doubleStain, "140616_doubleStainInvasion.txt", sep='\t')

# Plot Data
invD7NormP = ggplot(data=matchedCd133, aes(x=clone, y=mean, fill=cd133status)) + 
    scale_fill_manual(values=c("yellow", "skyblue3")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to control") +
    ggtitle("Invasive ability of GIC clones \nby CD133 status") +  # Set title
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


invD7DoubleStain = ggplot(data=doubleStain[!doubleStain$clone %in% "028",], aes(x=clone, y=mean, fill=cd133status)) + 
    #scale_fill_manual(values=c("firebrick", "skyblue3")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.9)) +
    xlab("Clone") + ylab("Surface area of gliomasphere \nnormalised to control") +
    ggtitle("Invasive ability of GIC clones \nby CD44 and CD133 status") +  # Set title
    theme_bw(base_size=16)#+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


multiplot(invD7StatsP, invD7NormP, invD7DoubleStain, cols=2)
