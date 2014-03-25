# This script will measure the proliferation over 3 days and also compare Resazurin and LDH
library(scatterplot3d)
source('~/Documents/Rscripts/qPCRFunctions.R')

setwd('~/Documents/Cell_biology/proliferation/Resazurin/140318_testing/linear/')

backgroundMeanSD <- function (dataFrame) {
    # Take the dataframe of raw data and then remove background fluorescence and take the mean and sd
    # Subtract the background (0 cells) which is the first row
    background = rowMeans(dataFrame[1,c(4:6)])
    # Subtract the background and bind back the metaData
    holder = dataFrame[,c(4:6)] - background
    result = cbind(dataFrame[,c(7:9)], holder)
    # Compute the mean and sd
    result$mean = rowMeans(result[,c(4:6)])
    result$sd = apply(result[,c(4:6)], 1, sd)
    return (result)
}

# The resazurin readings
day1 = read.delim('140319_day1_linearReplicates.txt')
day2 = read.delim('140320_day2_linearReplicates.txt')
day3 = read.delim('140321_day3_linearReplicates.txt')
resazurin = list(day1, day2, day3)

day1$day = 'one'
day2$day = 'two'
day3$day = 'three'

# the ldh readings
ldh5 = read.delim('../LDH/5min.txt')
ldh30 = read.delim('../LDH/30min.txt')

# Plot the replication correlations
par(mfrow=c(2,1))
s3d = scatterplot3d(day1$rep1, day1$rep2, day1$rep3, main="Day 1 Resazurin 3D Scatterplot", highlight.3d=TRUE, type="h", pch=16)
fit <- lm(rep1 ~ rep2+rep3, day1)
s3d$plane3d(fit)
t3d = scatterplot3d(day3$rep1, day3$rep2, day3$rep3, main="Day 3 Resazurin 3D Scatterplot", highlight.3d=TRUE, type="h", pch=16)
fit <- lm(rep1 ~ rep2+rep3, day3)
t3d$plane3d(fit)

# Remove background and take mean
day1Polish = backgroundMeanSD(day1)
day2Polish = backgroundMeanSD(day2)
day3Polish = backgroundMeanSD(day3)

# Make a big dataFrame of the Resazurin data
dataDaily = rbind(day1Polish[c(1:6),], day2Polish[c(1:6),], day3Polish[c(1:6),])
#write.table(dataDaily, '140325_resazurinSummary.txt', sep='\t')

# Plot Resazurin signal as a function of cell input
dailyPlot = ggplot(data=dataDaily, aes(x=cells, y=mean, group=day, colour=day)) + geom_line() + geom_point() +
                ggtitle('Cell input dependence on flourescence intensity') + xlab('Cells plated') + ylab('Fluorescence at 585nm') +
                theme_bw(base_size=16)
