setwd('~/Documents/eQTL/Matrix_eQTL_R/')

mutations = read.delim('130822_mutFileTestGoodNames.txt', row.names=1)
mutPatient = as.factor(colnames(mutations))
mutPatient = intersect(colnames(mutations), colnames(genes))
genes = read.delim('130822_geneFileTest.txt', row.names=1)

overlapGenes = genes[,mutPatient]
colnames(genes)
colnames(mutations)
colnames(overlapGenes)