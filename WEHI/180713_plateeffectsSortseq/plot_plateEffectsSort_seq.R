library(ggplot2)
library(reshape)
library(plyr)
library(platetools)

source("~/Code/Rscripts/Templates/multiplot.R")
mypath = "~/Documents/Lab_projects/Sort_seq/Data/Plots/"
setwd(mypath)

df = read.csv("longTableCountsERCC.csv")

first96 = df[c(1:1920),]
second96 = df[c(1921:3840),]
third96 = df[c(3841:5760),]
fourth96 = df[c(5761:7660),]

#### Make a dot plot vertical stacking of multiple groups
p1 <- ggplot(first96, aes(x=well, y=ercc)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))

p2 <- ggplot(second96, aes(x=well, y=ercc)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") + theme(legend.position="none") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))

p3 <- ggplot(third96, aes(x=well, y=ercc)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") + theme(legend.position="none") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))

p4 <- ggplot(fourth96, aes(x=well, y=ercc)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") + theme(legend.position="none") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))

multiplot(p1,p2,p3,p4,cols = 1)