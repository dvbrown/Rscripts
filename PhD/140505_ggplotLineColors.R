setwd('~/Dropbox/Demonstrating/2014/week7/')
library(ggplot2)
data = read.delim('humanAbPlot.txt', header=T)

dataPlot = ggplot(data=data, aes(x=dilution, y=mean, colour=Location)) + geom_line(aes(group=Location), size=1) + geom_point(size=3) + 
    scale_color_manual(values=c("darkgreen", "darkblue", 'firebrick4', 'orange2')) +
    ggtitle('Human blood samples by location discovered') + xlab('Dilution factor') + ylab('Mean isotype-adjusted absorbance') + 
    scale_y_continuous(breaks = round(seq(0, 1.75, by = 0.25), 2)) +
    theme_bw(base_size=20) + theme(plot.title = element_text(size = rel(1.25)))

dataPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1))