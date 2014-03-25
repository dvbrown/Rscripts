# This script will measure the proliferation over 3 days and also compare Resazurin and LDH
library(scatterplot3d)
source('~/Documents/Rscripts/qPCRFunctions.R')

setwd('~/Documents/Cell_biology/proliferation/Resazurin/140318_testing/linear/')
# The resazurin readings
day1 = read.delim('140319_day1_linearReplicates.txt')
day2 = read.delim('140320_day2_linearReplicates.txt')
day3 = read.delim('140321_day3_linearReplicates.txt')
resazurin = list(day1, day2, day3)

day1$cells = as.factor(day1$cells)
day2$cells = as.factor(day2$cells)
day3$cells = as.factor(day3$cells)

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