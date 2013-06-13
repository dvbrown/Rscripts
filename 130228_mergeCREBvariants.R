#A script to overlap a gene list with a melted vcf file (ie txt file)

variants = read.delim('~/Documents/RNAdata/RNAseqAnalysis/121121_TopHatNovelAlignment/varScanAnalysis/130211_reAnnotate/clone035allSamples.filterVarScan.fixVcf.annotateVariants.annotateDbSNP.snpSiftFilter.snpSiftExtract.txt')
genelist = read.delim('~/Documents/CREB/publicDataSets/Lembereger2008/GSE8948_RAW/GPL8321/results/130228_lembergerNormalvsKOcutLFC1.txt')

var = variants[,c(1,2,3,4,5,10,11,18)]

overlap = merge(var, genelist, by.x='EFF.0..GENE', by.y='GeneSymbol')

write.table(overlap,'~/Documents/CREB/130228_CREBvarScanOverlap.txt', sep='\t')