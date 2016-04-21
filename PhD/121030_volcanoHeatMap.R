library(DESeq)
library(ggplot2)
library(Heatplus)

finalResult = read.delim('121030_correctSamples/CD133lowCD133highCutoffs/121030_NaOmit.txt')

sigGenes=read.delim('121030_correctSamples/CD133lowCD133highCutoffs/121030_LFC1-5_Cutoff.txt', row.names=2)

require(ggplot2)
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
finalResult$threshold = as.factor(abs(finalResult$log2FoldChange) > 2 & finalResult$padj < 0.05)

##Construct the volcano plot object. First get the final result by running annotate ENSEMBL IDs script with genelist
g = ggplot(data=finalResult, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  opts(legend.position = "none", title=("Differential expression CD133 low vs CD133 high glioma initiating cells")
  ) +
    xlim(c(-5, 5)) + ylim(c(0, 50)) +
    xlab("log2 fold change") + ylab("-log10 FDR adjusted p-value")
g
#subset gene names for only significant genes
dd_text = finalResult[(abs(finalResult$log2FoldChange) > 4.2) & (finalResult$padj < 0.05),]
#add text to volcano
g + geom_text(data = dd_text, aes(x=log2FoldChange, y=-log10(padj),
                                  label=external_gene_id, size=0.25), colour="black")

#Draw a heatmap. First normalise the raw counts by library size
data = read.delim('121030_correctSamples/CD133lowCD133highCutoffs/121030_LFC1-5_Cutoff.txt')
row.names(data) = data$id
rawData = read.delim('121030_readCountsGenomeRanges.txt')
sigGenes = as.character(data$id) #use to index the raw data

dataMatrix = data.frame(row.names=colnames(rawData)
                        ,condition=c('CD133low_Rep1','CD133low_Rep2','CD133low_Rep3','CD133high_Rep1','CD133high_Rep1', 'CD133high_Rep1',
                                     'CD133unsort_rep1', 'CD133unsort_rep2', 'CD133unsort_rep3', 'test1','test2'),
                        libType=c('paired-end','paired-end','paired-end','paired-end','paired-end','paired-end',
                                  'paired-end','paired-end','paired-end','paired-end','paired-end'))
conds = dataMatrix$condition
countDataSet = newCountDataSet(rawData, conds)
#Get the library size scaled data
cds = estimateSizeFactors(countDataSet)
normaliseData = counts(cds, normalized=T) #use normalised data for the heat map
sigData = normaliseData[sigGenes,]
sigData = sigData[,c(1:9)]
geneNames = data$external_gene_id
row.names(sigData) = geneNames
#Now call the heat map
corrdist = function(x) as.dist(1-cor(t(x)))
heatmap(sigData, margins=c(7,5),cexRow=0.2)
