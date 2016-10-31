library(ggplot2)
library(RColorBrewer)
source('~/Code/Rscripts/Templates/multiplot.R')

inputDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/Counts_Maria/"
exportDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/"
setwd(inputDir)

# Read in samplesheet
samplesheet = read.csv("sampleSheetAll.csv", row.names = 1)
samplesheet = samplesheet[,c(1:4,6:8,9)]
rownames(samplesheet) = samplesheet[,5]
dat = read.delim("161031_librarySizeMaria.txt", row.names = 1)

batch1 = dat[dat$Batch == "A",]
platelay = data.frame (rown = rep (letters[1:8], 12), coln = rep (1:12, each = 8))
platelay1 = cbind(platelay, batch1)

ggplot(platelay1, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
  geom_point(aes(colour = GenesDetected, size= LibSize))  + theme_bw() + scale_size(range = c(5, 20)) +
  labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red")
