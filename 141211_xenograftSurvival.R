# Survival analysis of the big and failed mouse experient
library(survival)
library(coin)
library(ggplot2)

setwd("~/Documents/Xenograft/")
theData = read.delim("141211_mouseSurvival.txt")

colors = c("darkblue","orange", "red", "darkgreen", "yellow")
names(colors) = levels(theData$Sample)

data.surv = Surv(theData$Survival, event=theData$Event)
sur.fit = survfit(data.surv ~ Sample, theData)

plot(sur.fit, main='TCGA GBM cohort all patients classified by subtype',ylab='Survival probability', xlab='survival (days)', 
     col=colors,
     cex=1.75, conf.int=F, lwd=1.33, cex.axis=1.5, cex.lab=1.5)

legend('topright', names(colors), title="Coexpression subtype",
       col=colors,
       lwd=1.33, cex=1.75, bty='n', xjust=0.5, yjust=0.5)