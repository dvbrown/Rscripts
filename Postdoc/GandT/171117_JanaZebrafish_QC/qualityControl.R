library(ggplot2)
library(reshape)
library(xtabs)
source("~/Code/R_templates/multiplot.R")
setwd("~/Documents/Lab/Experiments/DNA_Quant_LibraryPrep/171109_JanaZebrafish/")
list.files()

dat = read.csv("171010_Jana_w_plateMap.csv",header=T,row.names = 1)
head(dat)
# Remove plate 4 which is empty
dat = dat[dat$Plate!= 4,]
dat = dat[!is.na(dat$Conc.x.Dil),]
dat = dat[!is.na(dat$Number),]
dat$Plate = as.factor(dat$Plate)
dat$DNAyield = dat$Conc.x.Dil * 30

### The number of samples that worked ###
dat$Worked = ifelse(dat$Conc.x.Dil > 0.9, TRUE, FALSE)

#### Plot boxplot of plate postion ####
byPlate = ggplot(dat, aes(x=factor(Plate), y=DNAyield, colour=Number)) + 
    geom_jitter() +
    xlab("Plate number") + ylab("DNA yield (ng)") + # Set axis labels
    ggtitle("DNA quantification Jana Zebrafish") + geom_hline(yintercept = 23, colour="red") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw(base_size=18) + theme(axis.text=element_text(size=18))
byPlate

#### Plot boxplot of sample ####
bySample = ggplot(dat, aes(x=factor(Identity), y=DNAyield, colour=Plate, shape=Number)) + 
   geom_jitter() +
    xlab("Sample") + ylab("DNA yield (ng)") + # Set axis labels
    ggtitle("DNA quantification Jana Zebrafish") + geom_hline(yintercept = 23, colour="red") + 
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(axis.text=element_text(size=18))
bySample