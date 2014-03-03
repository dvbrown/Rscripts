# This script will take the RNA-seq summary table and predict survival from CD133 status
source('~/Documents/Rscripts/140220_TCGAsignatureAnalysisFunctions.R')

setwd('~/Documents/RNAdata/')
data = read.delim('RNAseqProgress.txt')
data = data[c(1:14),]
data = data[c(1:11,13,14),]

correlation = lm(Survival_.month. ~ CD133pos_percent + Age, data=data)
#correlation = glm(Survival_.month. ~ CD133pos_percent + Age, data=data)

summary(correlation)

# Plot the scatter
termplot(correlation, terms='CD133pos_percent')
# Add fit lines

plot(data$CD133pos_percent, data$Survival_.month., main='Predictive value of CD133+ \npopulation with survival',
     xlab='Percent of population CD133+', ylab='Survival (months)', cex=1.33, col='black', pch=21, bg='blue')
abline(lm(data$Survival_.month. ~ data$CD133pos_percent), col="red") # regression line (y~x)

# Plot confidence intervals on line
# a <- predict(correlation, interval="confidence")
# lines(data$Survival_.month., a[,2], lty=2)
# lines(data$Survival_.month., a[,3], lty=2)

text(locator(1),labels='R squared = 0.04\np = 0.8', cex=0.75) #add the p-value to the graph