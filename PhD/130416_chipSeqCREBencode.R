setwd('~/Documents/CREB/ChIPseqENCODE/')
library(chipseq)
library(ShortRead)
library(BayesPeak)
A549_1 = read.delim('./wgEncodeHaibTfbsA549Creb1sc240V0416102Dex100nmPkRep1.broadPeak', header=F)

# function to import a MACS output file into a GRanges object
macs2GRange <-function(peaks) {
    myrange=GRanges(peaks[,1],IRanges(peaks[,2],peaks[,3], names=paste(peaks[,1],(peaks[,5]+peaks[,2]),sep=":")),
                    count=peaks[,6], score=peaks[,7],FE=peaks[,8],fdr=peaks[,9],summit=peaks[,5]
    )
    return(myrange)
}

x = macs2GRange(A549_1)