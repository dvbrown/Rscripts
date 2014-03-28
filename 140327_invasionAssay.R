# This script will measure the invasion of clone 035 across cells plated and the days observed
library(scatterplot3d)
source('~/Documents/Rscripts/qPCRFunctions.R')
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

meanSD <- function (dataFrame) {
    # Compute the mean and sd
    dataFrame$mean = rowMeans(dataFrame[,c(4:6)], na.rm=T)
    dataFrame$sd = apply(dataFrame[,c(4:6)], 1, sd, na.rm=T)
    return (dataFrame)
}

setwd('~/Documents/Cell_biology/microscopy/invasion/140222_clone035Test/linearReplicates/')
data = read.delim('140327_035_invasionResults.txt')

data$rep3 = as.numeric(data$rep3)

# Plot the replication correlations
par(mfrow=c(2,1))
s3d = scatterplot3d(data[c(5:8),4], data[c(5:8),5], data[c(5:8),6], main="Day 2 invasion assay replicates", highlight.3d=TRUE, type="h", pch=16,
                    xlab='replicate 1', ylab='replicate 2', zlab='replicate 3')
fit <- lm(rep1 ~ rep2+rep3, data[c(5:8),c(4:6)])
s3d$plane3d(fit)
t3d = scatterplot3d(data[c(13:16),4], data[c(13:16),5], data[c(13:16),6], main="Day 6 invasion assay replicates", highlight.3d=TRUE, type="h", pch=16,
                    xlab='replicate 1', ylab='replicate 2', zlab='replicate 3')
fit <- lm(rep1 ~ rep2+rep3, data[c(13:16),c(4:6)])
t3d$plane3d(fit)

# Take the mean and sd. Then remove the control
dataAv = meanSD(data)
dataMatrix = dataAv[dataAv$matrix %in% TRUE,]

# Plot Resazurin signal as a function of cell input
cellPlot = ggplot(data=dataMatrix, aes(x=cells, y=mean, group=day, colour=day)) + geom_line() + geom_point() +
    ggtitle('Area dependence on cell input') + xlab('Cells plated') + ylab('Surface area (uM)') + scale_fill_hue(name="Day plated") +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
    theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))
cellPlot

dataDaily$cells = as.factor(dataDaily$cells)
dailyPlot = ggplot(data=dataMatrix, aes(x=day, y=mean, group=cells, colour=cells)) + geom_line() + geom_point() +
    ggtitle('Invasion over time') + xlab('Day post invasion matrix') + ylab('Surface area (uM)') +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
    theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))
dailyPlot

multiplot(cellPlot, dailyPlot)

# Now compare no matrix to matrix
matrices = dataAv[dataAv$cells %in% 1000,]
matrixPlot = ggplot(data=matrices, aes(x=day, y=mean, group=matrix, colour=matrix)) + geom_line() + geom_point() +
    ggtitle('Invasion over time') + xlab('Day post invasion matrix') + ylab('Surface area (uM)') +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
    theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))
matrixPlot

matrices$day = as.factor(matrices$day)
bars = ggplot(data=matrices, aes(x=day, y=mean, fill=matrix)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_manual(values=c("skyblue", "pink")) +      # Set colors
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9)) +
    xlab("Days post invasion matrix addition") + ylab("Surface area (uM)") + # Set axis labels
    ggtitle("Invasion over time") +  # Set title
    theme_bw(base_size=20)
bars

multiplot(matrixPlot, bars)

write.table(dataAv, './analysis/140327_dataSummary.txt', sep='\t')