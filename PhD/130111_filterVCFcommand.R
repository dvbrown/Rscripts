args <- commandArgs(trailingOnly = TRUE)
#first argument is $input second $output is vcf file to filter and third argument is $output
rawVcf = read.delim(args[1],skip=115, header=F)
colnames(rawVcf) = c('#CHROM', 'POS','ID', 'REF','ALT','QUAL','FILTER', 'INFO', 'FORMAT', 'CD133n', 'sCD133n')
rawVcf$key = paste(rawVcf[,1], rawVcf$POS , sep="_")

cosmic = read.delim(args[2], skip=2,header=T)
cosmic$key = paste(cosmic[,1], cosmic$POS , sep="_")

filterVcf = merge.data.frame(cosmic, rawVcf, by.x='key', by.y='key')
filterVcf = filterVcf[,c(2:9, 18, 19,20)]
colnames(filterVcf) = c('#CHROM', 'POS','ID', 'REF','ALT','QUAL','FILTER', 'INFO', 'FORMAT', 'CD133n', 'sCD133n')

write.table(filterVcf, args[3] )