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
rawPlot = ggplot(data=day3Growth[day3Growth$treatment %in% growth,, 
                   aes(x=clone, y=ddCt, fill=cd133)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="sample") +      # Set legend title
    #coord_cartesian(ylim = c(0, 7.5)) +
    scale_y_continuous(breaks = round(seq(min(bindData$ddCt), max(bindData$ddCt), by = 1),1)) + # This modifies the scale of the y axis.
    xlab("Sample") + ylab("Gene expression normalised to CD133") + # Set axis labels
    ggtitle("Comapring CD133 status") +  # Set title
    theme_bw(base_size=20)
####################################################################################################################################



############################################## Read in the invasion assay readings #################################################
setwd('~/Documents/Cell_biology/microscopy/invasion/140414_invasion/')
invD3 = read.delim('140422_outputDay3.rep.txt')
invD7 = read.delim('140422_outputDay7.rep.txt')
####################################################################################################################################