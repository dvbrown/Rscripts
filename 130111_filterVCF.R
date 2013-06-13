source('/Users/d.brown6/Documents/Rscripts/120704-sortDataFrame.R')
rawVcf = read.delim('/Users/d.brown6/Documents/RNAdata/RNAseqAnalysis/121121_TopHatNovelAlignment/varScanAnalysis/clone035allSamples.sortSam.reorderSam.mpileup.pileup.varScan.vcf',skip=7, header=T)
#colnames(rawVcf) = c('#CHROM', 'POS','ID', 'REF','ALT','QUAL','FILTER', 'INFO', 'FORMAT', 'CD133n', 'sCD133n')
rawVcf$key = as.factor(paste(rawVcf[,1], rawVcf$POS , sep="_"))

cosmic = read.delim('/Users/d.brown6/Documents/public-datasets/annotationFiles/cosmic_vcfFile/CosmicMutantExport_v62_291112.tsv.txt',header=T)
cosmic1 = as.factor(gsub('-*', '', cosmic$Mutation.GRCh37.genome.position ))
cosmic$Mutation.GRCh37.genome.position = as.factor(cosmic1)
#cosmic = sort.dataframe(cosmic, cosmic$Mutation.GRCh37.genome.position, highFirst=F)

filterVcf = merge.data.frame(rawVcf, cosmic, by.x='key', by.y='Mutation.GRCh37.genome.position')
#filterVcf = merge.data.frame(rawVcf, cosmic, by=intersect(names(rawVcf$key), cosmic$Mutation.GRCh37.genome.position))

x = as.factor(intersect(rawVcf$key, cosmic$Mutation.GRCh37.genome.position))