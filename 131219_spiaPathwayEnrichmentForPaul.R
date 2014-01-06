# Measuring pathway enrichment using SPIA in RNA-seq batch1
library(SPIA)
library(pathview)
library(RColorBrewer)
library(ggplot2)
source('~/Documents/Rscripts/131218_ensemblToEnterezConversion.R')
source('~/Documents/Rscripts/120704-sortDataFrame.R')
setwd('~/Documents/CREB/paulEnrichment/')

zTransform = function(matrixElement, rowMean, rowSD ) { 
  # Convert to the unitless z-score based on a normal distribution
  z = (matrixElement - rowMean)/rowSD
  return (z)
}

############################# First z transform the TCGA data that Paul downloaded #############################

# Import the full dataset that Paul gave me. This is microarray data
data = read.delim('Colated TCGA_GBM.txt', row.names=1)

control = read.delim('colated TCGA_Controls.txt', row.names=1)
control = apply(control, 2, as.numeric)
# Remove the log tranformation if the data is still buggered
# control = log2(control)

dataNum = apply(data, 2, as.numeric)
# Remove the log tranformation if the data is still buggered
# dataNum = log2(dataNum)
head(dataNum)

# Obtain the z score for each gene
rowMean = rowMeans(control, na.rm=T)
rowStdDev = apply(control, 1, sd, na.rm=T)
dataNum = cbind(dataNum, rowMean)
# Compute the z-scores for the dataFrame
zScore = apply(dataNum, 2, zTransform, rowMean, rowStdDev)
row.names(zScore) = row.names(data)

# Don't log transform and see if that makes some sense
d = cbind(data, rowMean)
d = sort.dataframe(d, 596)

write.table(zScore, './131223_zTransormedTCGAgenes.txt', sep='\t', row.names=FALSE)

geneMean = rowMeans(zScore)
names(geneMean) = row.names(zScore)
zScore = cbind(zScore, geneMean)
x = as.data.frame(zScore)
x = sort.dataframe(x, 596, highFirst=T)

#Now build a the heat map
#geneMean = sort.default(geneMean)
#sub = geneMean[1:1000]
#zScorePlot = zScore[names(sub),]

#corrdist = function(x) as.dist(1-cor(t(x)))
#cc = brewer.pal(9, 'YlOrRd')
#heatmap(zScorePlot, col=cc, margins=c(7,5),cexRow=0.2, main='GBM gene expression TCGA', 
#        xlab='Patients', ylab='zTransformed genes')

####################################### Sumarise data and merge with CRE sites ########################################
# First try to subset untransformed data
ab = subset.data.frame(d, rowMean > 2)# , select=TCGA.02.0074.01A.01R.0195.07:TCGA.76.6280.01A.21R.1847.07)
ac = subset.data.frame(d, rowMean < -2)# , select=TCGA.02.0074.01A.01R.0195.07:TCGA.76.6280.01A.21R.1847.07)
abnormalGenes = rbind(ab, ac)
rm(ab, ac)

allIDs = read.delim('ensemblGeneIDsmart_export.txt')
row.names(allIDs) = allIDs$Ensembl.Gene.ID

# Extract the enterez IDs of all the genes tested
ensemblEnterezMap = ensembl2enterezConvert(allIDs)

# Merge the IDs with the untransformed data
abnormalGenes$gene = row.names(abnormalGenes)
mData = merge(abnormalGenes, ensemblEnterezMap, by.x='gene', by.y='SYMBOL')
#Just take the row mean
sumData = mData[,c(1, 597,598, 599)]
rm(d, allIDs, abnormalGenes, ensemblEnterezMap)

# Read in the CREB genes
crebTargets = read.delim('AllgenesCREBhits.txt')
sumData = merge.data.frame(sumData, crebTargets, by.x='ENTREZID', by.y='gene')

# Extract the interesting genes in GBM and CREB sites
intData = sumData[sumData$CREBcounts > 0,]

########################################## Now the SPIA enrichment ########################################

# Set up data for spia
genes <- abs(intData$rowMean) > 2
names(genes) <- intData$ENTREZID
all.genes <- crebTargets$gene

result.spia = spia(de=genes, all=all.genes, organism='hsa', nB=2000, plots=F)
result.spia[1:20, -12]

# tA is the equivaent of fold hange pertubation in the pathway.
# pXXX is all the various FDR corrections

plotP(result.spia, threshold=0.05)# x.lab='Enrichment score', ylab='Pertubation score', main='Disrupted pathways in short-term surviving GICs')

# Make a plot more amenable to publication
subsetResults.spia = result.spia[, c(1,2,3,4,6,7,9)]
colnames(subsetResults.spia) = c('Name', 'KEGG ID', 'Pathway size', 'Altered in GBM & CREB target genes', 'Pertubation_score', 'Pertubation score significance', 'Global_FDR_adjusted_pvalue')
subsetResults.spia[10,1] = ''
plot.default(subsetResults.spia$Pertubation_score, -log10(subsetResults.spia$Global_FDR_adjusted_pvalue),
             ylim=c(0,5))

################################# Construct the volcano plot object. ########################################

##Highlight genes that have an Pertubation score > 2 and a p-value < Bonferroni cut-off
subsetResults.spia$threshold = as.factor(abs(subsetResults.spia$Pertubation_score) > 2 & subsetResults.spia$Global_FDR_adjusted_pvalue < 0.05)

g = ggplot(data=subsetResults.spia, aes(x=Pertubation_score, y=-log10(Global_FDR_adjusted_pvalue), colour=threshold)) +
  geom_point(alpha=0.75, size=2) +
  opts(legend.position = "none", title=("Signaling Pathway Impact Analysis")
  ) +
  #xlim(c(-5, 5)) + ylim(c(0, 50)) +
  xlab("Pertubation score") + ylab("Significance")
g
#subset gene names for only significant genes
dd_text = subsetResults.spia[(abs(subsetResults.spia$Pertubation_score) > 2) & (subsetResults.spia$Global_FDR_adjusted_pvalue < 0.01),]
#add text to volcano
g + geom_text(data = dd_text, aes(x=Pertubation_score, y=-log10(Global_FDR_adjusted_pvalue),
                                  label=Name, size=0.2), colour="black")

########################################################################################################################

# output the results of the analysis
write.table(result.spia, './output/spia/140103_spiaResults.txt', sep='\t', row.names=F)
#write.table(spiaNoNAs, './output/spia/140103_spiaNoNAsResults.txt', sep='\t', row.names=F)

# Visualise the KEGG pathway data in R
# Use the results of plotP to get the pathway ids
setwd('./output/spia/')
pv.out = pathview(gene.data=intData$rowMean, pathway.id='04110', species='hsa', out.suffix='Cell cycle', kegg.native=T)
