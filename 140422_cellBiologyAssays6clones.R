source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data and then remove background fluorescence and take the mean and sd
    # Subtract the background (0 cells) which is the first row
    background = rowMeans(dataFrame[4,4:6], na.rm=T)
    # Subtract the background and bind back the metaData
    holder = dataFrame[,c(4:6)] - background
    result = cbind(dataFrame[,c(7:9)], holder)
    # Compute the mean and sd
    result$mean = rowMeans(result[,c(4:6)], na.rm=T)
    result$sd = apply(result[,c(4:6)], 1, sd, na.rm=T)
    return (result)
}

############################################## Read in the resazurin assay readings ###############################################
setwd('~/Documents/Cell_biology/proliferation/Resazurin/140417_6clones/analysis/')
growthD3 = read.delim('140414_day3_linearRep.txt')
growthD7 = read.delim('140417_day7_linearRep.txt')

par(mfrow=c(2,1))
plot(growthD7$rep1, growthD7$rep2, ylab='replicate 2', xlab='replicate1', main='consistency day7', pch=16)
abline(lm(growthD7$rep1~growthD7$rep2), col='red')
plot(growthD3$rep1, growthD3$rep2, ylab='replicate 2', xlab='replicate1', main='consistency day3', pch=16)
abline(lm(growthD3$rep1~growthD3$rep2), col='red')

# Subtract background and take mean and SD
day3Growth = backgroundMeanSD(growthD3)
day7Growth = backgroundMeanSD(growthD7)

# Plot the raw results
growthPlot3 = ggplot(data=day3Growth[day3Growth$treatment %in% 'growth',], 
                   aes(x=clone, y=mean, fill=cd133)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing proliferation at day 3 by CD133 status") +  # Set title
    theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))


growthPlot7 = ggplot(data=day7Growth[day7Growth$treatment %in% 'growth',], 
                     aes(x=clone, y=mean, fill=cd133)) + 
    scale_fill_manual(values=c("darkorange", "royalblue")) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
    xlab("Clone") + ylab("Fluorescent intensity") +
    ggtitle("Comparing proliferation at day 7 by CD133 status") +  # Set title
    theme_bw(base_size=20) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(growthPlot3, growthPlot7)
####################################################################################################################################



############################################## Read in the invasion assay readings #################################################
setwd('~/Documents/Cell_biology/microscopy/invasion/140414_invasion/')
invD3 = read.delim('140422_outputDay3.rep.txt')
invD7 = read.delim('140422_outputDay7.rep.txt')
####################################################################################################################################