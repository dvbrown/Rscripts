# TSS profile
library(GenomicAlignments)
library(ChIPseeker)
library(plyr)
library(Rsubread)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
setwd("~/Data/")

# Read in the bams
bams <- BamFileList(
  list.files(
    path = ("~/Data/"), 
    pattern = "^.*\\.bam$",
    full.names = TRUE))
bam = bams[[1]]

# Convert to genome ranges 
gas <- readGAlignments("atac_a3.bam")
grs <- granges(gas)

#promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
p <- promoterRegions(annotation = "hg38", upstream = 3000L, downstream = 3000L)
p <- as(p, "GRanges")
q = resize(p, fix = "center", 6001) # What does this mean?

tagMatrix <- getTagMatrix(grs, windows=q)
# tagHeatmap(tagMatrix, xlim=c(-1000, 3000), color="red")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

# Annotate the peaks to see if they fall in exons----
peakAnno <- annotatePeak(grs, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotDistToTSS(peakAnno, title="Distribution of reads loci relative to TSS")

# See if there is an overlap with the exons in the atac bam files ----
ex = makeGRangesFromDataFrame(getInBuiltAnnotation(annotation = "hg38"))
exResize = resize(ex, fix = "center", 4001) # What does this mean?
tagMatrix <- getTagMatrix(grs, windows=exResize)
plotAvgProf(tagMatrix, xlim=c(-2000, 2000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
