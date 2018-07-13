library(ggplot2)
library(reshape)
library(plyr)
library(platetools)

source("~/Code/Rscripts/Templates/multiplot.R")
mypath = "~/Documents/Lab_projects/Sort_seq/Data/180713_plateEffectsSORT_seq/"
setwd(mypath)
files = list.files(mypath)
sampleNames = substring(files, 1, 4)

genes <- as.character(read.table(files[1], header=TRUE, sep=",")[,1])

head(read.csv(files[1]))

df = do.call(cbind,lapply(files,function(fn)read.csv(fn)[,4]))

colnames(df) = sampleNames
df = as.data.frame(df)
df_log = df + 1
df_log = log2(df_log)

write.table(df, "summary.csv", sep=",", quote = F)
write.table(df_log, "../Plots/summaryLog2.csv", sep=",", quote = F)

# M19 is missing
well = num_to_well(1:384, plate = 384)
well = well[!well %in% "M19"]
df_log$well = well

long = as.data.frame(melt(df_log, id.vars = "well"))
long <- long[order(long$well),] 

# subset for testing

first96 = long[c(1:1920),]
second96 = long[c(1921:3840),]
third96 = long[c(3841:5760),]
fourth96 = long[c(5761:7660),]

#### Make a dot plot vertical stacking of multiple groups
p1 <- ggplot(first96, aes(x=well, y=value)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))
#ggsave("./Plots/boxfirst96.png", width=800, height=180, limitsize = F, units="mm")

p2 <- ggplot(second96, aes(x=well, y=value)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))
#ggsave("./Plots/boxsecond96.png", width=800, height=180, limitsize = F, units="mm")

p3 <- ggplot(third96, aes(x=well, y=value)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))
#ggsave("./Plots/boxthird96.png", width=800, height=180, limitsize = F, units="mm")

p4 <- ggplot(fourth96, aes(x=well, y=value)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(0.1)) + geom_jitter(aes(colour = variable), width = 0.1) +
  ggtitle("") +
  xlab("well") + ylab("Log2 Reads mapped to exons") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(text = element_text(size=16))
#ggsave("./Plots/boxfourth96.png", width=800, height=180, limitsize = F, units="mm")

multiplot(p1,p2,p3,p4,cols = 1)
#ggsave("./Plots/boxfirst96.png", width=800, height=180, limitsize = F, units="mm")
