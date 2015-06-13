library(ggplot2)
library(reshape)
library(plyr)
setwd("150612_hyperHypoxia/")

dat_all = read.csv("150612_hyperHypoxia.csv",row.names=1)
colnames(dat_all) = c("PDGC", "Treatment", "Viable", "CD44+/ CD133-",
                  "CD44+/ CD133+", "CD44-/ CD133+", "CD44-/ CD133+")

dat = dat_all[,c(1,2,4:7)]

# Convert from wide to long
mDat = melt(dat, id.vars = c("PDGC", "Treatment"))

# Biological replicates
bioRep = ddply(diff, .(PDGC, variable, Subpopulation), summarise, meanDiff = mean(value, na.rm=T), 
               sdDiff = sd(value, na.rm=T), reps=length(value))
bioRep$seDiff = bioRep$sdDiff / (bioRep$reps)


### Boxplot ###
box = ggplot(data=dat, aes(x=PDGC, y=ddCT)) +
    geom_boxplot() + geom_point(aes(colour=PDGC), size =3, alpha=0.7,  position = position_jitter(w = 0.175)) +
    ggtitle("qPCR Summary") + geom_hline(yintercept=0, colour="red") +
    xlab("Gene") + ylab("ddCt relative to CD133 negative") +
    scale_y_continuous(breaks = round(seq(-6, 6, by = 1),1)) +
    theme_bw(base_size=14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(text = element_text(size=20))
box