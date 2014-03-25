# This script will measure the proliferation over 3 days and also compare Resazurin and LDH
library(scatterplot3d)
source('~/Documents/Rscripts/qPCRFunctions.R')
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

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

day1$day_plated = '1'
day2$day_plated = '2'
day3$day_plated = '3'

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

########################################### Analysis of the daily measurements of Resazurin ##################################################################################### 
# Remove background and take mean
day1Polish = backgroundMeanSD(day1)
day2Polish = backgroundMeanSD(day2)
day3Polish = backgroundMeanSD(day3)

# Make a big dataFrame of the Resazurin data
dataDaily = rbind(day1Polish[c(1:6),], day2Polish[c(1:6),], day3Polish[c(1:6),])
#write.table(dataDaily, '140325_resazurinSummary.txt', sep='\t')

# Plot Resazurin signal as a function of cell input
cellPlot = ggplot(data=dataDaily, aes(x=cells, y=mean, group=day_plated, colour=day_plated)) + geom_line() + geom_point() +
                ggtitle('Cell input dependence on flourescence intensity') + xlab('Cells plated') + ylab('Fluorescence at 585nm') + scale_fill_hue(name="Day plated") +
                geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
                theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))
cellPlot

dataDaily$cells = as.factor(dataDaily$cells)
dailyPlot = ggplot(data=dataDaily, aes(x=day_plated, y=mean, group=cells, colour=cells)) + geom_line() + geom_point() +
    ggtitle('Time dependence on flourescence intensity') + xlab('Day of reading') + ylab('Fluorescence at 585nm') +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) +
    theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))
dailyPlot

multiplot(cellPlot, dailyPlot)
################################################################################################################################# 

########################################### Analysis temozolomide ###############################################################
tmz = day3[c(14:18),]
tmzSub = tmz[,4:6] - 1948
tmzSub$mean = rowMeans(tmzSub[,c(1:3)])
tmzSub$sd = apply(tmzSub[,c(1:3)], 1, sd)
tmzSub$conc = c(0, 6.25, 12.5, 25, 50)
tmz$conc = c(0, 6.25, 12.5, 25, 50)
tmz$mean = rowMeans(tmz[,c(4:6)])
tmz$sd = apply(tmz[,c(4:6)], 1, sd)

tmzPlot = ggplot(data=tmzSub, aes(x=conc, y=mean, group=1)) + geom_line(color='blue') + geom_point(color='darkblue') + 
            ggtitle('Effect of Temozolomide on cell number') + xlab('TMZ dose (uM)') + ylab('Fluorescence at 585nm') + 
            geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) + 
            theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))
tmzPlot

tmzRawPlot = ggplot(data=tmz, aes(x=conc, y=mean, group=1)) + geom_line(color='red') + geom_point(color='darkred') + 
    ggtitle('Effect of Temozolomide on cell number') + xlab('TMZ dose (uM)') + ylab('Fluorescence at 585nm') + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0)) + 
    theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))
tmzRawPlot
multiplot(tmzPlot, tmzRawPlot)
################################################################################################################################# 
########################################### Compare LDH to resazurin ###############################################################
ldh30$mean = rowMeans(ldh30[,c(4:6)])
resazurin = dataDaily[c(13:18),7]
ldh = ldh30$mean * 10000
# scale the resazurin and ldh. WORK ON THIS
pdf('../analysis/140325_correlationResazurinLDH.pdf', paper='a4')
plot(resazurin, ldh, main='Correlation Resazurin vs LDH')
abline(lm(resazurin~ldh), col="red") # regression line (y~x)
dev.off()
lines(lowess(resazurin, ldh), col="blue") # lowess line (x,y) 