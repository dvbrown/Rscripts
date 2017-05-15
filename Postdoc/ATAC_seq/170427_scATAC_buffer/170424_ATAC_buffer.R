library(ggplot2)
library(reshape)
library(plyr)

multmerge1 = function(mypath){
    filenames=list.files(path=mypath, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.delim(file=x,header=T, sep='\t', skip=6, nrows = 2)})
    return(datalist) }

multmerge2 = function(mypath){
    filenames=list.files(path=mypath, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.delim(file=x,header=F, sep='\t')})
    return(datalist) }

setwd('~/Data/ATAC/170420_ATACbuffer/LibMetrics/')
file_names = list.files()
source('~/Code/Rscripts/Templates/multiplot.R')

# Mung the library metrics
dat = multmerge1('~/Data/ATAC/170420_ATACbuffer/LibMetrics/')
df <- do.call("rbind", dat)
row.names(df) = list.files("./")

# Mung the mtDNA percentage
dat2 = multmerge2('~/Data/ATAC/170420_ATACbuffer/ChrsALignment/')
df2 <- do.call("cbind", dat2)
df3 = df2[,seq(2,73,2)]
row.names(df3) = df2[,1]
chrC = as.matrix(df3)
# Compute the sum
readSum = colSums(chrC)
percent = as.data.frame(chrC['chrM',] / readSum)
row.names(percent) = row.names(df)
colnames(percent) = "Percent_mtDNA"

# Bind the 2 together
datFrame = cbind(df, percent)
datFrame = datFrame[,c(2,4,5,7,10)]
rm(df, df2, df3, percent, chrC)
write.csv(datFrame, "~/Data/ATAC/170420_ATACbuffer/170426_atacMetrics.csv")
datFrame = read.csv("~/Data/ATAC/170420_ATACbuffer/170502_atacMetrics_annotate.csv", row.names = 1)
datFrame = subset(datFrame, Unmapped > 0)
datFrame$ReadsLog = log(datFrame$Unpaired_reads, 2)
datFrame$Reads10s = datFrame$Unpaired_reads / 1000

#### boxplot ####
color = rainbow(6)

p1 <- ggplot(datFrame, aes(factor(Lysis_Buffer), Percent_mtDNA)) +
    geom_boxplot() + geom_jitter(aes(colour=Cell_Line)) +
    ggtitle("Percent mtDNA ATAC-seq reads") + xlab("Lysis buffer") + ylab("Proportion mtDNA reads") +
    scale_colour_manual(values=c("blue", "red"), name="Lysis agent") +
    theme_bw(base_size=24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p1

p2 <- ggplot(datFrame, aes(factor(Lysis_Buffer), Reads10s)) +
    geom_boxplot() + geom_jitter(aes(colour=Cell_Line)) +
    ggtitle("") + xlab("Lysis buffer") + ylab("Unique reads (000s)") +
    scale_colour_manual(values=c("blue", "red"), name="Lysis agent") +
    theme_bw(base_size=24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p2

p3 <- ggplot(datFrame, aes(factor(Lysis_Buffer), Percent_unmapped)) +
    geom_boxplot() + geom_jitter(aes(colour=Cell_Line)) +
    ggtitle("") + xlab("Lysis buffer") + ylab("Percent unmapped reads") +
    scale_colour_manual(values=c("blue", "red"), name="Lysis agent") +
    theme_bw(base_size=24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p3

p3 <- ggplot(datFrame, aes(factor(Assay), Reads10s, fill = factor(Lysis_Buffer), label = factor(Percent))) +
    geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge") +
    scale_fill_manual(values=c("blue", "red"), name="Lysis buffer",
                      breaks=c("IGEPAL", "Digitonin"), labels=c("IGEPAL", "Digitonin")) +
    ggtitle("") + geom_text(check_overlap = TRUE, size = 4) +
    xlab("Lysis buffer") + ylab("Number of reads") + 
    theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=16))
p3

multiplot(p1,p2, cols=2)

#### Plot correlation of the replicate ####
correlation = ggplot(data=datFrame, aes(x=Reads10s, y=Percent_mtDNA, color=Lysis_Buffer)) + 
    geom_point(shape=19, size=4) + #geom_smooth(method=lm, colour='red') +
    xlab("Reads (1000s)") + ylab("Proportion mtDNA reads") + # Set axis labels
    ggtitle("") +  # Set title
    theme_bw(base_size=24)
correlation

#### Split samples by cell line andlysis buffer ####
GM12878 = datFrame[datFrame$Cell_Line == "GM12878",]
HCC38 = datFrame[!datFrame$Cell_Line == "GM12878",]

GM12878_Digit = GM12878[GM12878$Lysis_Buffer == "Digitonin",]
GM12878_IGEPAL = GM12878[!GM12878$Lysis_Buffer == "Digitonin",]
HCC38_Digit = HCC38[HCC38$Lysis_Buffer == "Digitonin",]
HCC38_IGEPAL = HCC38[!HCC38$Lysis_Buffer == "Digitonin",]


p5 <- ggplot(GM12878, aes(factor(Lysis_Buffer), Percent_mtDNA)) +
    geom_boxplot() + geom_jitter(aes(colour=Cell_Line)) +
    ggtitle("Percent mtDNA ATAC-seq reads") + xlab("Lysis buffer") + ylab("Proportion mtDNA reads") +
    scale_colour_manual(values=c("blue", "red"), name="Lysis agent") +
    theme_bw(base_size=24) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20))
p5