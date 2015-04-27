# This script will take the RNA-seq summary table and predict survival from CD133 status
source('~/Documents/Rscripts/140220_TCGAsignatureAnalysisFunctions.R')
library(ggplot2)

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

# Make a ggplot scatter
ggplot(data=data, aes(x=CD133pos_percent, y=Survival_.month., colour=Group)) + 
    scale_colour_manual(values=c("darkblue", "red")) +
    geom_point(shape=19, alpha=2/3, size=4) + geom_smooth(method=lm, colour='red', se=TRUE, alpha=1/9) +
    xlab("Percent of population expressing CD133") + ylab("Survival (months)") +
    ggtitle("") +  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    theme(text = element_text(size=20))

# Plot confidence intervals on line
# a <- predict(correlation, interval="confidence")
# lines(data$Survival_.month., a[,2], lty=2)
# lines(data$Survival_.month., a[,3], lty=2)
 
# text(locator(1),labels='R squared = 0.04\np = 0.8', cex=0.75) #add the p-value to the graph