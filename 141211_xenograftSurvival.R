# Survival analysis of the big and failed mouse experient
library(survival)
library(coin)
library(ggplot2)

setwd("~/Documents/Xenograft/")
theData = read.delim("141211_mouseSurvival.txt")

theDate = Sys.Date()

colors = c("darkblue","orange", "red", "darkgreen", "yellow")
names(colors) = levels(theData$Sample)

data.surv = Surv(theData$Survival, event=theData$Event)
sur.fit = survfit(data.surv ~ Sample, theData)

setwd("~/Documents/Xenograft/results/")
pdf(file=paste(theDate, "kmPlot.pdf", sep="_"), paper="a4r", useDingbats=F)

plot(sur.fit, main='Failed xenograft experiment',ylab='Survival probability', xlab='survival (days)', 
     col=colors,
     cex=1.75, conf.int=F, lwd=1.33, cex.axis=1.5, cex.lab=1.5)

legend('bottomleft', names(colors), title="Coexpression subtype",
       col=colors,
       lwd=1.33, cex=1.25, bty='n', xjust=0.5, yjust=0.5)
dev.off()

p = ggplot(theData, aes(Sample, Survival, fill=Sample)) +
    geom_boxplot() + geom_jitter() + scale_fill_manual(values=colors) +
    xlab("PDGC") + ylab("Survival (days)") +
    ggtitle("Failed xenograft experiment") + theme_bw(base_size=16) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24))


setwd("~/Documents/Xenograft/results/")
pdf(file=paste(theDate, "boxPlot.pdf", sep="_"), paper="a4r", useDingbats=F)
p
dev.off()

