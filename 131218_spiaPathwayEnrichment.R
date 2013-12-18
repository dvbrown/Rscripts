# Measuring pathway enrichment using SPIA in RNA-seq batch1
library(SPIA)
library(pathview)
source('~/Documents/Rscripts/131218_ensemblToEnterezConversion.R')
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/')

# Import the full dataset
data = read.delim('131021_shortVSlong.txt', row.names=2)
data = data[,c(2:8)]
data$ensembl = row.names(data)

# Add the enterez IDs
ensemblEnterezMap = ensembl2enterezConvert(data)
rnaseq.data.ent <- cbind(entid = ensemblEnterezMap$ENTREZID, data[ensemblEnterezMap$ENSEMBL,])

# Extract the differentially expressed genes
de.genes <- rnaseq.data.ent$logFC[data$FDR < 0.1]
names(de.genes) <- rnaseq.data.ent$entid[data$FDR < 0.1]
all.genes <- rnaseq.data.ent$entid

# Run spia
result.spia = spia(de=de.genes, all=all.genes, organism='hsa', nB=2000, plots=F)
result.spia$Name <- substr(spia.results$Name, 1, 25)
result.spia[1:20, -12]

# tA is the equivaent of fold hange pertubation in the pathway.
# pXXX is all the various FDR corrections

spiaNoNAs = result.spia[!is.na(result.spia$pPERT),]
plotP(spiaNoNAs, threshold=0.1)# x.lab='Enrichment score', ylab='Pertubation score', main='Disrupted pathways in short-term surviving GICs')

# output the results of the analysis
write.table(result.spia, './output/spia/131218_spiaResults.txt', sep='\t', row.names=F)
write.table(spiaNoNAs, './output/spia/131218_spiaNoNAsResults.txt', sep='\t', row.names=F)

# Visualise the KEGG pathway data in R
# Use the results of plotP to get the pathway ids
setwd('./output/spia/')
pv.out = pathview(gene.data=de.genes, pathway.id='04970', species='hsa', out.suffix='salivarySecretion', kegg.native=T)
pv.out = pathview(gene.data=de.genes, pathway.id='04512', species='hsa', out.suffix='ECMreceptor', kegg.native=T)
pv.out = pathview(gene.data=de.genes, pathway.id='05222', species='hsa', out.suffix='smallCellLungCancer', kegg.native=T)
pv.out = pathview(gene.data=de.genes, pathway.id='04020', species='hsa', out.suffix='calciumSignalling', kegg.native=T)
pv.out = pathview(gene.data=de.genes, pathway.id='04080', species='hsa', out.suffix='neuroactiveLigand', kegg.native=T)

