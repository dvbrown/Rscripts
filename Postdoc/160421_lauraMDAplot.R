library(ggplot2)
library(reshape)
library(plyr)
setwd("~/Documents/Presentations/2016/160421_ThierryRequest/")

dat = read.delim("lauraMDA.txt")

# Line plot
plt = ggplot(dat[dat$Sample %in% "Standard",], aes(x=Concentration_log10, y=Minutes, group=Sample)) +
  geom_point(shape=19, size = 3) + geom_smooth(method=lm, colour='darkblue', se=F) +
  ggtitle("") +  
  xlab("Concentration of DNA standard (log10)") + ylab("Time to threshold (minutes)") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=24))

plt + geom_point(data = dat[!dat$Sample %in% "Standard",], colour = "red", shape= 19, size = 4)
