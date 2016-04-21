#This makes a barplot with error bars
library(ggplot2)

coolPlot = function(data, sample, mean, stdDev) {
  newPlot = ggplot(data, aes(x=sample, y=mean)) +
                geom_bar(position=position_dodge(width=1), fill='blue', colour='black') +
                geom_errorbar(aes(ymin=mean-stdDev, ymax=mean+stdDev), 
                    width=.2,
                    position=position_dodge(.9), colour='black') +
                      ylab('Mean Ct') + opts(axis.text.x = theme_text(angle=90, hjust=1.2, size=12), 
                                             title='Threshold of primers that amplified')
  return (newPlot)
}