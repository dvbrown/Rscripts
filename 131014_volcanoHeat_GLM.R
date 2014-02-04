library(edgeR)
library(ggplot2)
library(Heatplus)
library(RColorBrewer)
setwd('~/Documents/RNAdata/danBatch1/bowtieGem/revHTSeq/GLMedgeR/140203_facsBatch/')

#finalResult = read.delim('GLMedgeR/131021_shortVSlong.txt')
finalResult = read.delim('140203_shortVSlong.txt')
#sigGenes=read.delim('GLMedgeR/131021_shortVSlongDEgenes.txt')
sigGenes=read.delim('140203_shortVSlongDEgenes.txt')

require(ggplot2)
##Highlight genes that have an absolute fold change > 1 and a p-value < Bonferroni cut-off
finalResult$threshold = as.factor(abs(finalResult$logFC) > 1 & finalResult$FDR < 0.1)

##Construct the volcano plot object. First get the final result by running annotate ENSEMBL IDs script with genelist
g = ggplot(data=finalResult, aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  opts(legend.position = "none", title=("Differential expression short vs long term survivors GIC RNA-seq batch1")
  ) +
  #xlim(c(-5, 5)) + ylim(c(0, 50)) +
  xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value")
g
#subset gene names for only significant genes
dd_text = finalResult[(abs(finalResult$logFC) > 1) & (finalResult$FDR < 0.00001),]
#add text to volcano
g + geom_text(data = dd_text, aes(x=logFC, y=-log10(FDR),
                                  label=external_gene_id, size=0.2), colour="black") +
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.title.y = element_text(size=14))


#Draw a heatmap. First normalise the raw counts by library size
data = read.delim('140203_shortVSlongDEgenes.txt')
row.names(data) = data$external_gene_id
rawData = read.delim('140203_normalisedCPM')
sigGenes = as.character(data$ensembl_gene_id) #use to index the raw data

sigData = as.matrix(rawData[sigGenes,])

#geneNames = data$external_gene_id ... FIX THIS

row.names(sigData) = geneNames
#Now call the heat map
corrdist = function(x) as.dist(1-cor(t(x)))
cc = brewer.pal(9, 'YlOrRd')
heatmap(sigData, col=cc, margins=c(7,5),cexRow=0.2, main='Hierarchical clustering GIC RNA-seq batch1', 
        xlab='Patient clones', ylab='Differentially expressed genes')
