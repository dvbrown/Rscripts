library(plyr)
library(ggplot2)
library(reshape)
source('~/Code/Rscripts/Templates/multiplot.R')

inputDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/Counts_Maria/"
exportDir = "/Users/u0107775/Data/GandT_Seq/161020_MariaXydia/"
setwd(inputDir)

# Read in samplesheet
samplesheet = read.csv("sampleSheetAll.csv", row.names = 1)
samplesheet = samplesheet[,c(1:4,6:8,9)]
rownames(samplesheet) = samplesheet[,5]

# Read in individaual batch GEMs
b1 = read.csv("Batch1/b1_htseq_counts_all.csv", row.names = 1)
b2 = read.csv("Batch2/b2_htseq_counts_all.csv", row.names = 1)
b3 = read.csv("Batch3/b3_htseq_counts_all.csv", row.names = 1)
b4 = read.csv("Batch4/b4_htseq_counts_all.csv", row.names = 1)

# Bind GEMS
count = cbind(b1, b2, b3, b4)
rm(b1, b2, b3, b4)
#write.table(count, "161026_geneExpMatrix.txt", sep="\t")

# Calculate library size and merge on sample sheet
libSize = as.data.frame(colSums(count))
row.names(libSize)
row.names(samplesheet)
# Calculate non-zero genes
geneSize =  as.data.frame(colSums(count != 0))
libSize = cbind(libSize, geneSize)
colnames(libSize) = c("LibSize", "NonZero")

libMatrix = merge.data.frame(samplesheet, libSize, by.x = 0, by.y = 0)

m = ggplot(libMatrix, aes(x=sample_ID, y=LibSize))  + geom_boxplot() +
  xlab("Sample") + ylab("Library Size") + # Set axis labels
  ggtitle("SMART-Seq lib size") + geom_jitter(aes(colour=Batch, width = 12),size=1.5)  +
  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size=20)) + guides(colour=FALSE)
m

n = ggplot(libMatrix, aes(x=sample_ID, y=NonZero))  + geom_boxplot() +
  xlab("Sample") + ylab("Genes with non-zero counts") + # Set axis labels
  ggtitle("SMART-Seq non-zero") + geom_jitter(aes(colour=Batch, width = 12),size=1.5)  +
  theme_bw(base_size=18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(text = element_text(size=20)) + guides(colour=FALSE)
n

multiplot(m, n, cols=2)