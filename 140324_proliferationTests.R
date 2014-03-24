# This script will measure the proliferation over 3 days and also compare Resazurin and LDH
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