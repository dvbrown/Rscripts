# Measuring pathway enrichment using SPIA in RNA-seq batch1
library(SPIA)
library(pathview)
source('~/Documents/Rscripts/131218_ensemblToEnterezConversion.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')
setwd('~/Documents/FredCSC/reformattedFiles/')

symbolEnterezMap = function(genesTestedforDE) {
    ensembl2enterez = select(org.Hs.eg.db, keys=(genesTestedforDE[,1]), cols=c('ENTREZID','ENSEMBL'), 
                             keytype='SYMBOL')
    ensembl2enterez = ensembl2enterez[!is.na(ensembl2enterez),]
    ensembl2enterez = ensembl2enterez[!is.na(ensembl2enterez$ENTREZID),]
    ensembl2enterez = ensembl2enterez[!duplicated(ensembl2enterez$ENTREZID),]
    return (ensembl2enterez)
}

# Import the full dataset
data = read.delim('130828_inputToSPIA.txt')
data = sort.dataframe(data, 5, highFirst=T)
# Remove duplicated entries should be 25,647 entries
dedupData = data[!duplicated(data[,1]),]
dupData = data[duplicated(data[,1]),]
dedupData$FDR = p.adjust(dedupData$P.value, method='BH')

# Add the enterez IDs
dataEnterez = symbolEnterezMap(dedupData)
data.ent = merge(dedupData, dataEnterez, by.x='Primary.Sequence.Name', by.y='SYMBOL')

# Extract the differentially expressed genes
de.genes = data.ent$log2FC[data.ent$FDR < 0.1]
names(de.genes) <- data.ent$ENTREZID[data.ent$FDR < 0.1]
all.genes <- data.ent$ENTREZID

# Run spia
result.spia = spia(de=de.genes, all=all.genes, organism='hsa', nB=2000, plots=F)
result.spia$Name <- substr(spia.results$Name, 1, 25)
result.spia[1:20, -12]

# tA is the equivaent of fold hange pertubation in the pathway.
# pXXX is all the various FDR corrections

spiaNoNAs = result.spia[!is.na(result.spia$pPERT),c(1:11)]
plotP(spiaNoNAs, threshold=0.1)#, x.lab='Enrichment score', ylab='Pertubation score', main='Disrupted pathways PXR overexpressing cells')
title(main='Pathway alteration in PXR \noverexpressing cells', xlab='Pathway enrichment', ylab='Pathway pertubation')

# output the results of the analysis
write.table(result.spia, './spia/140211_spiaResults.txt', sep='\t', row.names=F)
write.table(spiaNoNAs, './spia/140211_spiaNoNAsResults.txt', sep='\t', row.names=F)

# Visualise the KEGG pathway data in R
# Use the results of plotP to get the pathway ids
setwd('spia/')
pv.out = pathview(gene.data=de.genes, pathway.id='04010', species='hsa', out.suffix='MAPK', kegg.native=T)
pv.out = pathview(gene.data=de.genes, pathway.id='04142', species='hsa', out.suffix='Lysosome', kegg.native=T)
