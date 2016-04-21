library(edgeR)
setwd('~/Documents/public-datasets/firehose/stddata__2013_06_23/GBM/20130623/gdac.broadinstitute.org_GBM.mRNAseq_Preprocess.Level_4.2013062300.0.0/')

rawCounts = read.delim('GBM.uncv2.mRNAseq_raw_counts.txt', row.names=1)
scaledCounts = read.delim('GBM.uncv2.mRNAseq_scaled_estimate.txt')
normalised = read.delim('GBM.uncv2.mRNAseq_RSEM_normalized_log2.txt')

row.names(rawCounts)
dge = DGEList(counts=rawCounts)
scaledEdgeR = calcNormFactors(dge)

edgeRnorm = log2(scaledEdgeR$counts)

#leave the normalistion problem here for now. Get matrix eQTL working