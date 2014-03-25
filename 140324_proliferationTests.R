# This script will measure the proliferation over 3 days and also compare Resazurin and LDH
library(scatterplot3d)
source('~/Documents/Rscripts/qPCRFunctions.R')

setwd('~/Documents/Cell_biology/proliferation/Resazurin/140318_testing/linear/')

# The resazurin readings
day1 = read.delim('140319_day1_linearReplicates.txt')
day2 = read.delim('140320_day2_linearReplicates.txt')
day3 = read.delim('140321_day3_linearReplicates.txt')
resazurin = list(day1, day2, day3)

# day1$cells = as.factor(day1$cells)
# day2$cells = as.factor(day2$cells)
# day3$cells = as.factor(day3$cells)
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

backgroundMeanSD <- function (dataFrame) {
# Take the dataframe of raw data and then remove background fluorescence and take the mean and sd
  # Subtract the background (0 cells)
  background = rowMeans(dataFrame[1,c(4:6)])
  # Bind back the metaData
  holder = dataFrame[,c(4:6)] - background
  result = cbind(dataFrame[,c(7:11)], holder)
  # Compute the mean and sd
  result$mean = rowMeans(result[,c(4:6)])
  result$sd = apply(result[,c(4:6)], 1, sd)
  return (result)
}

# Subtract the background (0 cells)
day1Background = rowMeans(day1[1,c(4:6)])
day2Background = rowMeans(day2[1,c(4:6)])
day3Background = rowMeans(day3[1,c(4:6)])

day1Sub = day1[,c(4:6)] - day1Background
day1Sub = cbind(day1[,c(7:11)], day1Sub)

# Compute the mean and sd
day1$mean = rowMeans(day1[,c(4:6)])
day1$sd = apply(day1[,c(4:6)], 1, sd)
day2$mean = rowMeans(day2[,c(4:6)])
day2$sd = apply(day2[,c(4:6)], 1, sd)
day3$mean = rowMeans(day3[,c(4:6)])
day3$sd = apply(day3[,c(4:6)], 1, sd)

# Make a big dataFrame
data = rbind(day1, day2, day3)
#write.table(data, '140325_resazurinSummary.txt', sep='\t')

# Plot Resazurin signal as a function of cell input
ggplot(data=data, aes(x=cells, y=mean, group=day, colour=day)) + geom_line() + geom_point()

day1 = ggplot(data=rawData[!rawData$gene.x %in% c('GFAP'),], aes(x=origin.x, y=ddCt_B2M_020_P, fill=gene.x)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") + 
    scale_fill_hue(name="Gene") +      # Set legend title
    xlab("Sample") + ylab("ddCt") + # Set axis labels
    ggtitle("B2M as a house keeping gene") +  # Set title
    theme_bw(base_size=14)