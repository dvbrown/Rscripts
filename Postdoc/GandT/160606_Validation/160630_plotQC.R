library(plyr)
library(ggplot2)
library(reshape)

inputDir = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/QC/"
exportDir = "/Users/u0107775/Data/GandT_Seq/160606_ValidationGandT/DNA/QC/Rplot/"

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

ggplot(mdataInteresting, aes(x=variable, y=value, label=LIBRARY)) + geom_boxplot() +
  #scale_colour_manual(values=variable) +
  xlab("Metric") + ylab("Number") + # Set axis labels
  ggtitle("G&T-Seq") + geom_point(aes(colour=variable)) +
  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size=20)) + guides(colour=FALSE)

ggplot(percentages, aes(x=variable, y=value, label=LIBRARY)) + geom_boxplot() +
  #scale_colour_manual(values=variable) +
  xlab("Metric") + ylab("Number") + # Set axis labels
  ggtitle("G&T-Seq") + geom_point(aes(colour=variable)) + geom_text() +
  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size=20)) + guides(colour=FALSE)

ggplot(mdata[mdata$variable %in% c("ESTIMATED_LIBRARY_SIZE"),], aes(x=variable, y=value, label=LIBRARY)) + geom_boxplot() +
  #scale_colour_manual(values=variable) +
  xlab("Metric") + ylab("Number") + # Set axis labels
  ggtitle("G&T-Seq") + geom_point(aes(colour=variable)) + geom_text() +
  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size=20)) + guides(colour=FALSE)
