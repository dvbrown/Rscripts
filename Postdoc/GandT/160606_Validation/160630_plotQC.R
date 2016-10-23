library(plyr)
library(ggplot2)
library(reshape)
source('~/Code/Rscripts/Templates/multiplot.R')

inputDir = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/QC/"
exportDir = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/QC/Rplot/"

# Get a list of the bedtools output files you'd like to read in
setwd(inputDir)
print(files <- list.files(pattern="*.bam.txt$"))

# Create lists to hold information for each alignment,
# and read the data into these lists.
qc <- list()
qc_cumul <- list()
for (i in 1:length(files)) {
  qc[[i]] <- read.delim(files[i], skip=6, header=T, nrows=1)
}

df <- ldply(qc, data.frame)
mdata <- melt(df, id=c("LIBRARY"))
mdataInteresting = mdata[!mdata$variable %in% c("ESTIMATED_LIBRARY_SIZE","READ_PAIR_OPTICAL_DUPLICATES", "UNPAIRED_READ_DUPLICATES", "PERCENT_DUPLICATION"),]
percentages = mdata[mdata$variable %in% c("PERCENT_DUPLICATION"),]

a = ggplot(mdataInteresting, aes(x=variable, y=value, label=LIBRARY)) + geom_boxplot() +
      #scale_colour_manual(values=variable) +
      xlab("") + ylab("Number") + # Set axis labels
      ggtitle("G&T-Seq") + geom_point(aes(colour=variable)) +
      theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      theme(text = element_text(size=20)) + guides(colour=FALSE)

d = ggplot(mdataInteresting, aes(x=variable, y=value, label=LIBRARY)) + geom_boxplot() +
  #scale_colour_manual(values=variable) +
  xlab("") + ylab("Number") + # Set axis labels
  ggtitle("G&T-Seq") + geom_point(aes(colour=variable)) +
  coord_cartesian(ylim = c(0, 4000000)) +
  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size=20)) + guides(colour=FALSE)

b = ggplot(percentages, aes(x=variable, y=value, label=LIBRARY)) + geom_boxplot() +
    #scale_colour_manual(values=variable) +
    xlab("") + ylab("Duplicate rate") + # Set axis labels
    ggtitle("G&T-Seq") + geom_point(aes(colour=variable)) + geom_text() +
    theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(text = element_text(size=20)) + guides(colour=FALSE)

c = ggplot(mdata[mdata$variable %in% c("ESTIMATED_LIBRARY_SIZE"),], aes(x=variable, y=value, label=LIBRARY)) + geom_boxplot() +
      #scale_colour_manual(values=variable) +
      xlab("") + ylab("Number") + # Set axis labels
      ggtitle("G&T-Seq") + geom_jitter(aes(colour=variable, width = 0.1),size=2)  +
      theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      theme(text = element_text(size=20)) + guides(colour=FALSE)

multiplot(a, d, b, c, cols=2)
write.table(mdata, "data.txt", sep='\t')

# Get the MAPD scores
mapd = read.delim("mapdScoreAnnotated.txt")
mapd = mapd[c(2:6)]

m = ggplot(mapd, aes(x=Cell_type, y=MAPD))  + geom_boxplot(colour=Cell_number) +
  xlab("") + ylab("MAPD DNA analysis") + # Set axis labels
  ggtitle("MAPD") + geom_jitter(aes(colour=Cell_number, width = 12),size=3,)  +
  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size=20)) + guides(colour=FALSE)
m