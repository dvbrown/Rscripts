library(ggplot2)
source('~/Code/Rscripts/Templates/multiplot.R')

plot96wellMap <- function (sampleSubset, variable_2_colour, variable_2_size) {
  # DATAFRAME    sampleSubset = the 96 samples you want to plot. This must be in order with well A12 = sample 12 and well B1 = sample 13
  # STRING     variable_2_colour = the measurement you want mapped to colour
  # STRING     variable_2_size = the measurement you want mapped to circle size
  
  wellMap = data.frame (rown = rep (LETTERS[1:8], each=12), coln = rep (1:12, 8))
  platelay = cbind(wellMap, sampleSubset)
  
  # ggplot using well coordinates generate dabove
  plt = ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
    geom_point(aes_string(colour = variable_2_colour, size= variable_2_size)) + scale_size(range = c(5, 10)) +
    labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red") +
    theme_bw() + theme(axis.text=element_text(size=16))
  return(plt)
}

b1 = plot96wellMap(batch1, "GenesDetected", "LibSize")
b1

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
batch2 = datMerge[c(97:192),]
batch3 = datMerge[c(193:288),]
batch4 = datMerge[c(289:384),]

# Insert well positions batch 1
wellMap = data.frame (rown = rep (LETTERS[1:8], each=12), coln = rep (1:12, 8))
platelay = cbind(wellMap, batch1)
# ggplot using well coordinates
b1 = ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
    geom_point(aes(colour = GenesDetected, size= LibSize)) + scale_size(range = c(5, 10)) +
    labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red") +
    theme_bw() + theme(axis.text=element_text(size=16))

# Insert well positions batch 2
platelay = cbind(wellMap, batch2)
# ggplot using well coordinates
b2 = ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
  geom_point(aes(colour = GenesDetected, size= LibSize)) + scale_size(range = c(5, 10)) +
  labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red") +
  theme_bw() + theme(axis.text=element_text(size=16))

# Insert well positions batch 3
platelay = cbind(wellMap, batch3)
# ggplot using well coordinates
b3 = ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
  geom_point(aes(colour = GenesDetected, size= LibSize)) + scale_size(range = c(5, 10)) +
  labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red") +
  theme_bw() + theme(axis.text=element_text(size=16))

# Insert well positions batch 4
platelay = cbind(wellMap, batch4)
# ggplot using well coordinates
b4 = ggplot(platelay, aes(y = factor(rown, rev(levels(rown))),x = factor(coln))) + 
  geom_point(aes(colour = GenesDetected, size= LibSize)) + scale_size(range = c(5, 10)) +
  labs(x=NULL, y = NULL) + scale_color_gradient(low="blue", high="red") +
  theme_bw() + theme(axis.text=element_text(size=16))

# PLOT ALL PLATES
multiplot(b1, b2, b3, b4, cols = 2)