# This script will measure the invasion of clone 035 across cells plated and the days observed
library(scatterplot3d)
source('~/Documents/Rscripts/qPCRFunctions.R')
source('~/Documents/Rscripts/140211_multiplotGgplot2.R')

setwd('~/Documents/Cell_biology/microscopy/invasion/140222_clone035Test/linearReplicates/')
data = read.delim('140327_035_invasionResults_edit.txt')

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