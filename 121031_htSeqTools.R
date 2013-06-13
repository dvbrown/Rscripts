#Compute some coverage statistics
library(htSeqTools)
library(rtracklayer)
setwd('~/Documents/RNAdata/RNAseqAnalysis/bedFiles/')
neg = import.bed('CD133n_B.bed') #read in the BED files. These were generated from aligned bam files.
pos = import.bed('CD133p_B.bed')
combine = RangedDataList(neg, pos) #A container for multiple bed files

cmds1 = cmds(combine, k=2) #generate a multi-dimensonal scaling plot to determine similarity. eg like PCA
plot(cmds1, main='MDS plot variability in CD133 sorted cells')
ssdCoverage(combine) #Measure similarity in coverage
giniCoverage(combine[['neg']], mk.plot=T, chrLengths='integer') #reports an error that bed fileshave null and missing values
filterDuplReads(combine)
