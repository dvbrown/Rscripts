#merge the logFC data with a bed file to make a nice heatmap
logFC = read.delim('~/Documents/RNAdata/RNAseqAnalysis/121105_trimmomaticReads/mergedBam/121107_mergeSortTopHatAlignIndex/121113_DEseqNoreplicates/130320_CD133lowVSCD133Ed.txt')

bed = read.delim('~/Documents/RNAdata/RNAseqAnalysis/121105_trimmomaticReads/mergedBam/121107_mergeSortTopHatAlignIndex/121113_DEseqNoreplicates/DEgenes.bed.txt')
be = unique(bed)

data = merge.data.frame(be, logFC, by.x='name2', by.y='external_gene_id')

write.table(data, '~/Bioinformatics/circos-0.63-4/dan/130320_DEgenesHeat.txt',sep='\t')