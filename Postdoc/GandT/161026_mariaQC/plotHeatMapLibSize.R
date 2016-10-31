library(ggplot2)
library(RColorBrewer)
source('~/Code/Rscripts/Templates/multiplot.R')

inputDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/Counts_Maria/"
exportDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/"
setwd(inputDir)

# Read in samplesheet
samplesheet = read.csv("sampleSheetAll.csv", row.names = 6)
samplesheet = samplesheet[,c(1:4,6:8,9)]
#rownames(samplesheet) = samplesheet[,5]
dat = read.delim("161031_librarySizeMaria.txt", row.names = 1)

# Sort samples by well position using sample sheet
datMerge = merge.data.frame(dat, samplesheet, by.x = 0, by.y = 0)
datMerge <- datMerge[order(datMerge$sample_number),]

# Get the first batch of samples sample_number 1 - 96
batch1 = datMerge[c(1:96),]

# Insert well positions
wellMap = data.frame (rown = rep (letters[1:8], each=12), coln = rep (1:12, 8))
platelay = cbind(wellMap, batch1)

# ggplot using well coordinates
ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
  geom_point(aes(colour = GenesDetected, size= LibSize))  + theme_bw() + scale_size(range = c(5, 20)) +
  labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red")
