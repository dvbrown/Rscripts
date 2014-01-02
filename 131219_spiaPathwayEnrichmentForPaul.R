# Measuring pathway enrichment using SPIA in RNA-seq batch1
library(SPIA)
library(pathview)
library(RColorBrewer)
source('~/Documents/Rscripts/131218_ensemblToEnterezConversion.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')
setwd('~/Documents/CREB/paulEnrichment/')

zTransform = function(matrixElement, rowMean, rowSD ) { 
  # Convert to the unitless z-score based on a normal distribution
  z = (matrixElement - rowMean)/rowSD
  return (z)
}

# LOG TRANSFORM

############################# First z transform the TCGA data that Paul downloaded #############################

# Import the full dataset that Paul gave me. This is microarray data
data = read.delim('Colated TCGA_GBM.txt', row.names=1)

control = read.delim('colated TCGA_Controls.txt', row.names=1)
control = apply(control, 2, as.numeric)
# Remove the log tranformation if the data is still buggered
control = log2(control)

dataNum = apply(data, 2, as.numeric)
# Remove the log tranformation if the data is still buggered
dataNum = log2(dataNum)
head(dataNum)

# Obtain the z score for each gene
rowMean = rowMeans(control, na.rm=T)
rowStdDev = apply(control, 1, sd, na.rm=T)
dataNum = cbind(dataNum, rowMean)
# Compute the z-scores for the dataFrame
zScore = apply(dataNum, 2, zTransform, rowMean, rowStdDev)
row.names(zScore) = row.names(data)
# Check the calculation of the z score
x$rowMean

write.table(zScore, './131223_zTransormedTCGAgenes.txt', sep='\t', row.names=FALSE)

geneMean = rowMeans(zScore)
names(geneMean) = row.names(zScore)
zScore = cbind(zScore, geneMean)
x = as.data.frame(zScore)
x = sort.dataframe(x, 596, highFirst=T)

#Now build a the heat map
geneMean = sort.default(geneMean)
sub = geneMean[1:1000]
zScorePlot = zScore[names(sub),]

corrdist = function(x) as.dist(1-cor(t(x)))
cc = brewer.pal(9, 'YlOrRd')
heatmap(zScorePlot, col=cc, margins=c(7,5),cexRow=0.2, main='GBM gene expression TCGA', 
        xlab='Patients', ylab='zTransformed genes')

############################# Now the SPIA enrichment #############################

# y = rowMedians(zScore[,c(1:596)], na.rm=T)
# abnormalGenes = x[(x$geneMean > 2),]
abnormalGenes = subset.data.frame(x, -2 < geneMean > 2, select=TCGA.02.0074.01A.01R.0195.07:TCGA.76.6280.01A.21R.1847.07)

allIDs = read.delim('ensemblGeneIDsmart_export.txt')
row.names(allIDs) = allIDs$Ensembl.Gene.ID

# Extract the enterez IDs of all the genes tested
ensemblEnterezMap = ensembl2enterezConvert(allIDs)

# Extract the interesting genes in GBM and CREB sites
de.genes <- data$CREBcounts[data$CREBcounts > 0]
names(de.genes) <- data$gene_id[data$CREBcounts > 0]
all.genes <- ensemblEnterezMap$ENTREZID

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

